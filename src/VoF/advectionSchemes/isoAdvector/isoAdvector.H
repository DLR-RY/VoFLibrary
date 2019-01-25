/*---------------------------------------------------------------------------*\
    Modified work | Copyright (c) 2017-2019, German Aerospace Center (DLR)
-------------------------------------------------------------------------------
License
    This file is part of the VoFLibrary source code library, which is an 
	unofficial extension to OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    isoAdvector

Description
    Calculates the new VOF (alpha) field after time step dt given the initial
    VOF field and a velocity field U and face fluxes phi. The fluid transport
    calculation is based on an idea of using isosurfaces to estimate the
    internal distribution of fluid in cells and advecting such isosurfaces
    across the mesh faces with the velocity field interpolated to the
    isosurfaces.

    Note: normal vector points in the fluid phase

    Reference:
    \verbatim
        Roenby, J., Bredmose, H. and Jasak, H. (2016).
        A computational method for sharp interface advection
        Royal Society Open Science, 3
        doi 10.1098/rsos.160405

        Henning Scheufler, Johan Roenby,
        Accurate and efficient surface reconstruction from volume
        fraction data on general meshes,
        Journal of Computational Physics, 2019,
        doi 10.1016/j.jcp.2019.01.009

    \endverbatim

    Original code supplied by Johan Roenby, DHI (2016)


SourceFiles
    isoAdvector.C

\*---------------------------------------------------------------------------*/

#ifndef isoAdvector_H
#define isoAdvector_H

#include "typeInfo.H"
#include "advectionSchemes.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"

#include "volPointInterpolation.H"
#include "reconstructionSchemes.H"
#include "cutFaceAdvect.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace advection
{

/*---------------------------------------------------------------------------*\
                Class isoAdvector Declaration
\*---------------------------------------------------------------------------*/

class isoAdvector
    : public advectionSchemes
{

private:
    // Private data types

        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Face volumetric transport
        surfaceScalarField dVf_;

        //- An isoCutFace object to get access to its face cutting functionality
        cutFaceAdvect advectFace_;

        // Switches and tolerances. Tolerances need to go into toleranceSwitches

            //- Number of alpha bounding steps
            label nAlphaBounds_;

            //- Tolerance for search of isoFace giving specified VOF value
            scalar vof2IsoTol_;

            //- Tolerance for marking of surface cells:
            //  Those with surfCellTol_ < alpha1 < 1 - surfCellTol_
            scalar surfCellTol_;
            
            //- Switch controlling whether to use isoface normals for interface
            //  orientation (default corresponding to false) to base it on 
            //  a smoothed gradient of alpha calculation (giving better results 
            //  on tri on tet meshes).
            bool gradAlphaBasedNormal_;
            
            //- Print isofaces in a <case>/isoFaces/isoFaces_#N.vtk files.
            //  Intended for debugging
            bool writeIsoFacesToFile_;

        // Cell and face cutting

            //- List of surface cells
            DynamicLabelList surfCells_;

            //- Bool list for cells that have been touched by the bounding step
            boolList cellIsBounded_;

            //- True for all surface cells and their neighbours
            boolList checkBounding_;

            //- Storage for boundary faces downwind to a surface cell
            DynamicLabelList bsFaces_;

            //- Storage for boundary surface iso face centre
            DynamicVectorList bsx0_;

            //- Storage for boundary surface iso face normal
            DynamicVectorList bsn0_;

            //- Storage for boundary surface iso face speed
            DynamicScalarList bsUn0_;

            //- Storage for boundary surface iso value
      //      DynamicScalarList bsf0_;

            //- Smallest face area
            scalar minMagSf_;

            scalar limiter_;


        // Additional data for parallel runs

            //- List of processor patch labels
            DynamicLabelList procPatchLabels_;

            //- For each patch if it is a processor patch this is a list of the
            //  face labels on this patch that are downwind to a surface cell.
            //  For non-processor patches the list will be empty.
            List<DynamicLabelList> surfaceCellFacesOnProcPatches_;

            //- For each face calculate volumetric face transport during dt
            void timeIntegratedFlux();

            //- Calculate volumetric face transport during dt given the isoFace
            //  data provided as input for face facei

            //- For a given cell return labels of faces fluxing out of this cell
            //  (based on sign of phi)
            void setDownwindFaces
            (
                const label ci,
                DynamicLabelList& downwindFaces
            ) const;

            // Limit fluxes
            void limitFluxes();

            void limitFluxes(const volScalarField::Internal& Sp,const volScalarField::Internal& Su);

            // Bound fluxes
            void boundFromAbove
            (
                const scalarField& alpha1,
                surfaceScalarField& dVfcorrected,
                DynamicLabelList& correctedFaces
            );

            // Bound fluxes
            void boundFromAbove
            (
                const scalarField& alpha1,
                surfaceScalarField& dVfcorrected,
                DynamicLabelList& correctedFaces,
                const volScalarField::Internal& Sp,
                const volScalarField::Internal& Su
            );

            // Bound fluxes
	    // is need as the source terms need to be modified if the inverse equation is solved: alpha2=(1-alpha)
	    // to avoid this alpha1 will be bounded to >= 0 not  alpha2 <= 1
            void boundFromBelow
            (
                const scalarField& alpha1,
                surfaceScalarField& dVfcorrected,
                DynamicLabelList& correctedFaces,
                const volScalarField::Internal& Sp,
                const volScalarField::Internal& Su
            );

            //- Given the face volume transport dVf calculates the total volume
            //  leaving a given cell. Note: cannot use dVf member because
            //  netFlux is called also for corrected dVf
            scalar netFlux
            (
                const surfaceScalarField& dVf,
                const label celli
            ) const;

            //- Determine if a cell is a surface cell
            bool isASurfaceCell(const label celli) const
            {
                return
                (
                    surfCellTol_ < alpha1_[celli]
                 && alpha1_[celli] < 1 - surfCellTol_
                );
            }

            //- Clear out isoFace data
            void clearIsoFaceData()
            {
                surfCells_.clear();
                bsFaces_.clear();
                bsx0_.clear();
                bsn0_.clear();
                bsUn0_.clear();
                if(mesh_.topoChanging())
                {
                    // Introduced resizing to cope with changing meshes
                    checkBounding_.resize(mesh_.nCells());
                    cellIsBounded_.resize(mesh_.nCells());
                }
                checkBounding_ = false;
                cellIsBounded_ = false;
            }


        // Face value functions needed for random face access where the face
        // can be either internal or boundary face

            //- Return face value for a given Geometric surface field
            template<typename Type>
            Type faceValue
            (
                const GeometricField<Type, fvsPatchField, surfaceMesh>& f,
                const label facei
            ) const;

            //- Set face value for a given Geometric surface field
            template<typename Type>
            void setFaceValue
            (
                GeometricField<Type, fvsPatchField, surfaceMesh>& f,
                const label facei,
                const Type& value
            ) const;


        // Parallel run handling functions

            //- Synchronize dVf across processor boundaries using upwind value
            void syncProcPatches
            (
                surfaceScalarField& dVf,
                const surfaceScalarField& phi
            );

            //- Check if the face is on processor patch and append it to the
            //  list of surface cell faces on processor patches
            void checkIfOnProcPatch(const label facei);

        void calcFaceFlux();

        //- Apply the bounding based on user inputs
        void applyBruteForceBounding();


private:

    // Private Member Functions

        //- Disallow copy construct
        isoAdvector(const isoAdvector&);

        //- Disallow default bitwise assignment
        void operator=(const isoAdvector&);



public:

    //- Runtime type information
    TypeName("isoAdvector");


    // Constructors

        //- Construct from components
        isoAdvector
        (
            volScalarField &alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U
        );


    //- Destructor
    virtual ~isoAdvector();


    //- advection of the interface
    virtual void advect();

    //- advection of the interface with source terms
    virtual void advect(const volScalarField::Internal& Sp,const volScalarField::Internal& Su);


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End advection Foam
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
