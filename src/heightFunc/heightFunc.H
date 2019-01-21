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
    Foam::heightFunc

Description
    Calculates a reconstructed distance function next to interface,
    the gradient in interface cells  and interpolates normal vector
    in interface cells

Author
    Henning Scheufler, DLR

SourceFiles
    heightFunc.C

\*---------------------------------------------------------------------------*/

#ifndef heightFunc_H
#define heightFunc_H

#include "OFstream.H"
#include "cutCellIso.H"
#include "fvMesh.H"
#include "globalIndex.H"
#include "scalarMatrices.H"
#include "surfaceFields.H"
#include "token.H"
#include "volFields.H"

#include "cellNeighbours.H"
#include "leastSquareCell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class heightFunc Declaration
\*---------------------------------------------------------------------------*/

class heightFunc
{


private:

    // Private data

    // reference to mesh
    const fvMesh& mesh_;

    // reference to volume fraction field
    const volScalarField& alpha1_;

    // value it which it is considered an interface
    scalar surfCellTol_;

    // is true if a cell has a interface
    boolList interface_;

    // is true of one adajecent point is next to an interface
    boolList nextToInterface_;

    // is true if point is next to an interface cell
    boolList interfacePoint_;

    // provides stencil and map
    cellNeighbours nei_;

    // checks if cell has interface
    bool isASurfaceCell(const scalar alpha1) const
    {
         return
         (
             surfCellTol_ < alpha1
             && alpha1 < (1 - surfCellTol_)
         );
    }

    // loop Stencil and calc weighted distance to interface
    void loopStencil
    (
        const DynamicList < label >& neiCells,
        const Map < List <vector> >& map,
        const globalIndex& globalIdx,
        const volVectorField& centre,
        const volVectorField& normal,
        const point& p,
        scalar& averageDist,
        scalar& avgWeight
    );

    // loop Stencil and find nearestInterface cell(nearestNeighbour interpolation)
    void loopStencilFindNearNeighbour
    (
        const DynamicList < label >& neiCells,
        const Map < List <vector> >& map,
        const globalIndex& globalIdx,
        const volVectorField& centre,
        const volVectorField& normal,
        const point& p,
        vector& normalNearNei,
        DynamicList<vector>& interfaceNormals
    );

    // loop Stencil and find nearestInterface cell(nearestNeighbour interpolation)
    void loopStencilInterpolateNormal
    (
        const DynamicList < label >& neiCells,
        const Map < List <vector> >& map,
        const globalIndex& globalIdx,
        const volVectorField& centre,
        const volVectorField& normal,
        const point& p,
        vector& normalNearNei,
        DynamicList<vector>& interfaceNormals
    );

    // loop Stencil and find nearestInterface cell(nearestNeighbour interpolation)
    void loopStencilEstimatedNormalAndCentre
    (
        const DynamicList < label >& neiCells,
        const Map < List <vector> >& map,
        const globalIndex& globalIdx,
        const volVectorField& centre,
        const volVectorField& normal,
        const point& p,
        vector& avgNormal,
        vector &interpolCentre
    );

    // loop Stencil and find nearestInterface cell(nearestNeighbour interpolation)
    void loopStencilAdvectHeightAndEstimateNormal
    (
        const DynamicList < label >& neiCells,
        const Map < List <vector> >& map,
        const globalIndex& globalIdx,
        const volVectorField& centre,
        const volVectorField& normal,
        const volVectorField& U,
        const point& p,
        scalar& height,
        vector& avgNormal,
        vector &interpolCentre
    );

public:

    // Constructors

    //- Construct from fvMesh and a scalarField
    //  Length of scalarField should equal number of mesh points
    heightFunc(const fvMesh&, const volScalarField& alpha1 ,const scalar& tol);


    // Member functions

    //- updates interfacePoint_ interface_ and  nextToInterface_
    void markCellsNearSurf();

    //- calculates height function
    void calcHeightFunc(const volVectorField& centre, const volVectorField& normal, volScalarField& height);

    //- interpolates normals
    void interpolateNormals(const volVectorField& centre, const volVectorField& normal, volVectorField& surfnormals);

    //- interpolates normals
    void interpolateNormals(const volVectorField& centre, const volVectorField& normal, const volVectorField& U, volVectorField& surfnormals);

    //- calculates the gradient
    void grad(const volScalarField& height, volVectorField& surfnormals);

    //- calculates the gradient
    void grad(const volScalarField& height,const volVectorField& centre, volVectorField& surfnormals);

    void calcResidual(const volVectorField& normal, const volVectorField& interfaceNormal, Map<scalar>& normalResidual, Map<scalar>& avgAngle);

    void calcResidual(const volVectorField& normal, Map<vector>& oldNormal, Map<scalar>& normalResidual);

    //- interpolate a point value
    void interpolatePoints(const volScalarField& height, scalarField& ap);

    //- return interfacePoints_
    const boolList& interfacePoints()
    {
        return interfacePoint_;
    }

    //- returns boolList of cells with an interface
    const boolList& interfaceCells()
    {
        return interface_;
    }

    //- returns boolList of cells with an interface and the neighbour cells
    const boolList& nextToInterface()
    {
        return nextToInterface_;
    }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
