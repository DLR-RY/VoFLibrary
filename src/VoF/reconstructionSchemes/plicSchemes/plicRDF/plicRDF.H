/*---------------------------------------------------------------------------*\
            Copyright (c) 2017-2019, German Aerospace Center (DLR)
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
    Foam::plicRDF

Description
    Reconstructs an interface (centre and normal vector) consisting of planes
    to match the internal fluid distribution in cells. The interface normals
    are estimated by least square gradient scheme on the RDF function (height).
    Uses the normal from the previous times step as intial guess.

    Reference:
    \verbatim
        Henning Scheufler, Johan Roenby,
        Accurate and efficient surface reconstruction from volume
        fraction data on general meshes,
        Journal of Computational Physics, 2019,
        doi 10.1016/j.jcp.2019.01.009

    \endverbatim

Author
    Henning Scheufler, DLR, all rights reserved.

SourceFiles
    plicRDF.C

\*---------------------------------------------------------------------------*/

#ifndef plicRDF_H
#define plicRDF_H

#include "typeInfo.H"
#include "reconstructionSchemes.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"


#include "surfaceIteratorPLIC.H"
#include "reconstructedDistanceFunction.H"
#include "zoneDistribute.H"




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reconstruction
{

/*---------------------------------------------------------------------------*\
                Class plicRDF Declaration
\*---------------------------------------------------------------------------*/

class plicRDF
    : public reconstructionSchemes
{

private:

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Interpolation object from cell centres to points
        //Map<vector> interfaceNormal_;
        DynamicField<vector> interfaceNormal_;

        



        
    // Switches and tolerances. Tolerances need to go into toleranceSwitches


        //- Tolerance for search of isoFace giving specified VOF value
        scalar vof2IsoTol_;

        //- Tolerance for marking of surface cells:
        //  Those with surfCellTol_ < alpha1 < 1 - surfCellTol_
        scalar surfCellTol_;

        //  tolerance
        scalar tol_;

        //  relative tolerance
        scalar relTol_;

        // number of times that the interface ist reconstructed has to be bigger than 2
        label iteration_;

        // interpolated normal from previous time step
        bool interpolateNormal_;


        // calculates the RDF function
        reconstructedDistanceFunction RDF_;

        // provides stencil and map
        zoneDistribute& exchangeFields_;

        //- surfaceIterator finds the plane centre for specified VOF value
        surfaceIteratorPLIC sIterPLIC_;


private:

    // Private Member Functions

        //- Disallow copy construct
        plicRDF(const plicRDF&);

        //- Disallow default bitwise assignment
        void operator=(const plicRDF&);

        void setInitNormals(bool interpolate);

        void gradSurf(const volScalarField& phi);

        void calcResidual
        (
            Map<scalar>& normalResidual,
            Map<scalar>& avgAngle
        );

        void interpolateNormal();



public:

    //- Runtime type information
    TypeName("plicRDF");


    // Constructors

        //- Construct from components
        plicRDF
        (
                volScalarField& alpha1,
                const surfaceScalarField& phi,
                const volVectorField& U,
                dictionary& dict
        );


    //- Destructor
    virtual ~plicRDF();

    // reconstruct  interface
    virtual void reconstruct();

    virtual void mapAlphaField();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End reconstructionSchemes Foam
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
