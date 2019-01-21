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
    Foam::leastSquareCell

Description
    Estimates the gradient with a least square scheme in a cell

Author
    Henning Scheufler, DLR

SourceFiles
    leastSquareCell.C

\*---------------------------------------------------------------------------*/

#ifndef leastSquareCell_H
#define leastSquareCell_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "OFstream.H"
#include "cutCellIso.H"
#include "scalarMatrices.H"
#include "globalIndex.H"

#include "simpleMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class leastSquareCell Declaration
\*---------------------------------------------------------------------------*/

class leastSquareCell
{



private:
    // Private data

        //- Mesh whose cells and faces to cut at their intersection with an
        //  isosurface.
    const fvMesh& mesh_;

    // checks if 2D
    //
    // -1 is 3D
    // 1 sym plane is normal to x
    // 2 sym plane is normal to y
    // 3 sym plane is normal to z
    label twoDim_;

    // reference to globalNumbering
    const globalIndex& globalNumbering_;

    //- gives patchNumber and patchFaceNumber for a given Geometric volume field
   void patchIFace
   (
       const label faceI,
       label& patchNumber,
       label& patchFaceNumber
   ) const;


    bool isInternalPoint(const label& pI) // maybe inline ?
    {
        return pI < mesh_.nInternalPoints();
    }

    // fill matrix values for 2D matrix
    void fill2DMatrix
    (
            const vector& vec,
            const scalar& height,
            simpleMatrix<scalar>& A

    );

    // fill matrix values for 3D matrix
    void fill3DMatrix
    (
            const vector& vec,
            const scalar& height,
            simpleMatrix<scalar>& A

    );

    // estimates the gradient in 2D
    vector grad2D
    (
            const label& cellI,
            const vector& pos,
            const volScalarField& phi,
            const DynamicList < label > & stencil,
            const Map < scalar >& map,
            const Map < vector >& mapCC
    );

    // estimates the gradient in 3D
    vector grad3D
    (
            const label& cellI,
            const vector& pos,
            const volScalarField& phi,
            const DynamicList < label > & stencil,
            const Map<scalar> &map,
            const Map < vector >& mapCC
    );

    // interpolate a point value
    scalar interpolate2D
    (
        const label& pI,
        const volScalarField& phi,
        const Map <List <scalar> >& map,
        const Map <List <vector> >& mapCC
    );


    // interpolate a point value
    scalar interpolate3D
    (
        const label& pI,
        const volScalarField& phi,
        const Map <List <scalar> >& map,
        const Map <List <vector> >& mapCC
    );


public:

    // Constructors

        //- Construct from fvMesh and globalIndex
        leastSquareCell(const fvMesh&, const globalIndex& gblIdx);


    // Member functions

        //- returns the gradient of the cell
        vector calcGrad
        (
                const label& cellI,
                const vector& pos,
                const volScalarField& phi,
                const DynamicList<label> &stencil,
                const Map <scalar> &map,
                const Map <vector>& mapCC
        );
       // vector calcGrad(const label& cellI, const volScalarField& phi, const volVectorField& centre); // interpolate a point value
        //- interpolate a point value
        scalar interpolate
        (
            const label& pI,
            const volScalarField& phi,
            const Map <List <scalar> >& map,
            const Map <List <vector> >& mapCC
        );



};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
