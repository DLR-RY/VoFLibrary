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


\*---------------------------------------------------------------------------*/

#include "leastSquareCell.H"

#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"
#include "wedgePolyPatch.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::leastSquareCell::patchIFace
(
    const label faceI,
    label& patchNumber,
    label& patchFaceNumber
) const
{
    if (mesh_.isInternalFace(faceI))
    {
        patchNumber = -1;
    }
    else
    {
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

        // Boundary face. Find out which face of which patch
        const label patchI = pbm.patchID()[faceI - mesh_.nInternalFaces()];

        if (patchI < 0 || patchI >= pbm.size())
        {
            FatalErrorInFunction << "Cannot find patch for face " << faceI
                                 << abort(FatalError);
        }

        // Handle empty patches

        const polyPatch& pp = pbm[patchI];
        if (isA<emptyPolyPatch>(pp) || pp.empty() || isA<wedgePolyPatch>(pp) ||
            isA<processorPolyPatch>(pp))
        {
            patchNumber = -1;
        }
        else
        {
            patchNumber = patchI;
        }

        const label patchFaceI = pp.whichFace(faceI);
        patchFaceNumber = patchFaceI;
    }
}

void Foam::leastSquareCell::fill2DMatrix
(
    const vector& vec,
    const scalar& height,
    simpleMatrix< scalar >& A
)
{
    scalar cpX = 0;
    scalar cpZ = 0;

    if (twoDim_ == 0) // very stupid solution
    {
        cpX = vec.y();
        cpZ = vec.z();
    }
    if (twoDim_ == 1)
    {
        cpX = vec.x();
        cpZ = vec.z();
    }
    if (twoDim_ == 2)
    {
        cpX = vec.x();
        cpZ = vec.y();
    }

    A[0][0] += cpX * cpX;
    A[0][1] += cpX * cpZ;
    A[0][2] += cpX;

    A[1][0] += cpX * cpZ;
    A[1][1] += cpZ * cpZ;
    A[1][2] += cpZ;

    A[2][0] += cpX;
    A[2][1] += cpZ;
    A[2][2] += 1;

    A.source()[0] += cpX * height;
    A.source()[1] += cpZ * height;
    A.source()[2] += height;
}

void Foam::leastSquareCell::fill3DMatrix
(
    const vector& vec,
    const scalar& height,
    simpleMatrix< scalar >& A
)
{
    A[0][0] += vec.x() * vec.x();
    A[0][1] += vec.x() * vec.y();
    A[0][2] += vec.x() * vec.z();
    A[0][3] += vec.x();

    A[1][0] += vec.y() * vec.x();
    A[1][1] += vec.y() * vec.y();
    A[1][2] += vec.y() * vec.z();
    A[1][3] += vec.y();

    A[2][0] += vec.z() * vec.x();
    A[2][1] += vec.z() * vec.y();
    A[2][2] += vec.z() * vec.z();
    A[2][3] += vec.z();

    A[3][0] += vec.x();
    A[3][1] += vec.y();
    A[3][2] += vec.z();
    A[3][3] += 1;

    A.source()[0] += vec.x() * height;
    A.source()[1] += vec.y() * height;
    A.source()[2] += vec.z() * height;
    A.source()[3] += height;
}

Foam::vector Foam::leastSquareCell::grad2D
(
    const label& cellI,
    const vector& pos,
    const volScalarField& phi,
    const DynamicList< label >& stencil,
    const Map< scalar >& map,
    const Map< vector >& mapCC
)
{

    simpleMatrix< scalar > A(3, 0, 0); // 3 mal 3 matrix
    scalarField coeffs(3);             // coeffs

    forAll(stencil, i)
    {
        const label gIdx = stencil[i];
        if (globalNumbering_.isLocal(gIdx))
        {
            // also includes boundary conditions
            const label localCellI = globalNumbering_.toLocal( gIdx);
            if (localCellI < mesh_.nCells())
            {

                vector vec = (mesh_.C()[localCellI] - pos);

//                if(pos != mesh_.C()[cellI] && i == 0)
//                {

//                    fill2DMatrix(vector::zero, 0, A);
//                }
//                else
//                {
                    fill2DMatrix(vec, phi[localCellI], A);
//                }
            }
//            else
//            {
//                // face is saved in map
//                const scalar height = map[gIdx];
//                const vector c = mapCC[gIdx];

//                vector vec = (c - pos);

//                fill2DMatrix(vec, height, A);
//            }
        }
        else
        {
            const scalar height = map[gIdx];
            const vector c = mapCC[gIdx];
            vector vec = (c - pos);

            fill2DMatrix(vec, height, A);
        }
    }

    coeffs = A.LUsolve();

    vector normal(0, 0, 0);
    if (twoDim_ == 0) // very stupid solution
    {

        normal.y() = coeffs[0];
        normal.z() = coeffs[1];
    }
    if (twoDim_ == 1)
    {
        normal.x() = coeffs[0];
        normal.z() = coeffs[1];
    }
    if (twoDim_ == 2)
    {
        normal.x() = coeffs[0];
        normal.y() = coeffs[1];
    }

    return normal;
}

Foam::vector Foam::leastSquareCell::grad3D
(
    const label& cellI,
    const vector& pos,
    const volScalarField& phi,
    const DynamicList< label >& stencil,
    const Map< scalar >& map,
    const Map< vector >& mapCC
)
{
    simpleMatrix< scalar > A(4, 0, 0); // 3 mal 3 matrix
    scalarField coeffs(4);             // coeffs

    //    fill3DMatrix(vector(0,0,0),phi[cellI],A);

    forAll(stencil, i)
    {
        const label gIdx = stencil[i];
        if (globalNumbering_.isLocal(gIdx))
        {
            // also includes boundary faces
            const label localCellI = globalNumbering_.toLocal(gIdx);
            if (localCellI <  mesh_.nCells())
            {
                vector vec = (mesh_.C()[localCellI] - pos);

                if(pos != mesh_.C()[cellI] && i == 0)
                {
                    fill3DMatrix(vector::zero, 0, A);
                }
                else
                {
                    fill3DMatrix(vec, phi[localCellI], A);
                }


            }
//            else
//            {
//                // cell is saved in map
//                //                if(!map.found(gIdx))
//                //                {
//                //                    Pout << "proc " <<
//                //                    globalNumbering_.whichProcID (gIdx) <<
//                //                    endl;
//                //                    map[gIdx];
//                //                }
//                const scalar height = map[gIdx];
//                const vector c = mapCC[gIdx];

//                vector vec = (c - pos);

//                fill3DMatrix(vec, height, A);
//            }
        }
        else
        {
            const scalar height = map[gIdx];
            const vector c = mapCC[gIdx];
            //            if(!map.found(gIdx))
            //            {
            //                Pout << "proc " << globalNumbering_.whichProcID
            //                (gIdx) << endl;
            //                map[gIdx];
            //            }

            vector vec = (c - pos);

            fill3DMatrix(vec, height, A);
        }
    }

    coeffs = A.LUsolve();

    vector normal(coeffs[0], coeffs[1], coeffs[2]);

    return normal;
}

Foam::scalar Foam::leastSquareCell::interpolate2D
(
    const label& pI,
    const volScalarField& phi,
    const Map <List <scalar> >& map,
    const Map <List <vector> >& mapCC
)
{
    simpleMatrix< scalar > A(3, 0, 0); // 3 mal 3 matrix
    scalarField coeffs(3);             // coeffs

    const labelListList& cPoints = mesh_.pointCells();


    if(map.found(pI))
    {
        forAll(map[pI],i)
        {
            const scalar height = map[pI][i];
            const vector c = mapCC[pI][i];

            vector vec = (c - mesh_.points()[pI]);

            fill2DMatrix(vec, height, A);
        }

    }
    else
    {
        forAll(cPoints[pI],j) // loop over all cells attached to the point
        {
            const label& pCellI = cPoints[pI][j];
           // f = a*x + b*z + c
            vector vec = mesh_.C()[pCellI] - mesh_.points()[pI];

            fill2DMatrix(vec, phi[pCellI], A);

        }
    }

    coeffs = A.LUsolve(); // gauss with piviting
    return coeffs[2];
}

Foam::scalar Foam::leastSquareCell::interpolate3D
(
    const label& pI,
    const volScalarField& phi,
    const Map <List <scalar> >& map,
    const Map <List <vector> >& mapCC
)
{
    simpleMatrix< scalar > A(4, 0, 0); // 3 mal 3 matrix
    scalarField coeffs(4);             // coeffs

    const labelListList& cPoints = mesh_.pointCells();
    const labelListList& fPoints = mesh_.pointFaces();
    if(map.found(pI))
    {

        forAll(map[pI],i)
        {
            const scalar height = map[pI][i];
            const vector c = mapCC[pI][i];

            vector vec = (c - mesh_.points()[pI]);

            fill3DMatrix(vec, height, A);
        }
    }
    else
    {
        forAll(cPoints[pI],j) // loop over all cells attached to the point
        {
            const label& pCellI = cPoints[pI][j];
           // f = a*x + b*z + c
            vector vec = mesh_.C()[pCellI] - mesh_.points()[pI];

            fill3DMatrix(vec, phi[pCellI], A);

        }

//        forAll(fPoints[pI],j) // loop over all cells attached to the point
//        {
//            const label& faceI = fPoints[pI][j];
//           // f = a*x + b*z + c
//            label patchi = -1;
//            label i = -1;

//            //faceValue<Type>(f, faceI, patchi, i);
//            patchIFace(faceI, patchi, i);
//            if (patchi != -1)
//            {
//                scalar val = phi.boundaryField()[patchi][i];
//                vector vec = mesh_.C().boundaryField()[patchi][i] - mesh_.points()[pI];
//                fill3DMatrix(vec, val, A);

//            }





//        }


    }

    coeffs = A.LUsolve(); // gauss with piviting
    return coeffs[3];

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::leastSquareCell::leastSquareCell
(
    const fvMesh& mesh,
    const globalIndex& gblIdx
)
:
    mesh_(mesh),
    twoDim_(-1),
    globalNumbering_(gblIdx)
{

    label dimensions = 0;
    for (direction cmpt = 0; cmpt < vector::nComponents; cmpt++)
    {
        dimensions += pos(mesh_.geometricD()[cmpt]) * mesh_.geometricD()[cmpt];
    }

    if (dimensions == 3)
    {
        twoDim_ = -1;
    }
    else if (dimensions == 2)
    {
        if (mesh_.geometricD().x() == -1)
        {
            twoDim_ = 0;
        }
        else if (mesh_.geometricD().y() == -1)
        {
            twoDim_ = 1;
        }
        else if (mesh_.geometricD().z() == -1)
        {
            twoDim_ = 2;
        }
    }
    else
    {
        // throw error
    }

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::leastSquareCell::calcGrad
(
    const label& cellI,
    const vector& pos,
    const volScalarField& phi,
    const DynamicList< label >& stencil,
    const Map< scalar >& map,
    const Map< vector >& mapCC
)
{
    if (twoDim_ != -1)
    {
        return grad2D(cellI, pos,phi, stencil, map, mapCC);
    }
    else
    {
        return grad3D(cellI, pos,phi, stencil, map, mapCC);
    }
}

Foam::scalar Foam::leastSquareCell::interpolate
(
    const label& pI,
    const volScalarField& phi,
    const Map <List <scalar> >& map,
    const Map <List <vector> >& mapCC
)
{



    if(twoDim_ != -1)
    {
        return interpolate2D(pI,phi,map,mapCC);
    }
    else
    {
        return interpolate3D(pI,phi,map,mapCC);
    }
}

// ************************************************************************* //
