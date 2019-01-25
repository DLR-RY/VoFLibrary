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
    Foam::surfaceIteratorPLIC

Description
    Finds the cutValue that matches the volume fraction



Author
    Johan Roenby, DHI, all rights reserved.
    Modified Henning Scheufler, DLR

SourceFiles
    surfaceIteratorPLIC.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceIteratorPLIC_H
#define surfaceIteratorPLIC_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "OFstream.H"
#include "cutCellPLIC.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class surfaceIteratorPLIC Declaration
\*---------------------------------------------------------------------------*/

class surfaceIteratorPLIC
{




private:
    // Private data

        //- Mesh whose cells and faces to cut at their intersection with an
        //  isosurface.
        const fvMesh& mesh_;

        // cuts the cell and returns volume centre and normal
        cutCellPLIC cutCell_;//(mesh_,normal);

        //- Tolerance for marking of surface cells:
        //  Those with surfCellTol_ < alpha1 < 1 - surfCellTol_
        scalar surfCellTol_;


public:

    // Constructors

        //- Construct from fvMesh and a scalarField
        //  Length of scalarField should equal number of mesh points
        surfaceIteratorPLIC(const fvMesh&,const scalar& tol);


    // Member functions
        //- Determine if a cell is a surface cell
        bool isASurfaceCell(const scalar alpha1) const
        {
            return
            (
                surfCellTol_ < alpha1
             && alpha1 < 1 - surfCellTol_
            );
        }






        //- finds matching cutValue for the given value fraction
        //- returns the cellStatus
        label vofCutCell
        (
            const label celli,
            const scalar alpha1,
            const scalar tol,
            const label maxIter,
            vector normal
        );

        //- returns the centre point of cutted volume
        point subCellCentre()
        {
            return cutCell_.subCellCentre();
        }


        //- returns the volume of cutted volume
        scalar subCellVolume()
        {
            return cutCell_.subCellVolume();
        }

        //- returns the centre of cutting isosurface
        point surfaceCentre()
        {
            return cutCell_.faceCentre();
        }

        //- returns the area vector of cutting isosurface
        vector surfaceArea()
        {
            return cutCell_.faceArea();
        }


        //- Volume of Fluid for cellI (subCellVolume_/mesh_.V()[cellI])
        scalar VolumeOfFluid()
        {
            return cutCell_.VolumeOfFluid();
        }

        //- returns the cutValue
        scalar cutValue() const
        {
            return cutCell_.cutValue();
        }

        //- returns the cellStatus
        label cellStatus()
        {
            return cutCell_.cellStatus();
        }

        //- returns the points of the cutting isosurface in sorted order
        DynamicList<point> facePoints()
        {
            return cutCell_.facePoints();
        }


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
