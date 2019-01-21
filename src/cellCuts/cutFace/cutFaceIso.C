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
    
\*---------------------------------------------------------------------------*/

#include "cutFaceIso.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cutFaceIso::cutFaceIso(const fvMesh& mesh, scalarField& f)
:
    cutFace(mesh),
    mesh_(mesh),
    f_(f),
    subFaceCentre_(vector::zero),
    subFaceArea_(vector::zero),
    subFacePoints_(10),
    surfacePoints_(4),
    pointStatus_(10),
    weight_(10),
    faceStatus_(-1)
{
    clearStorage();
}

// * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::cutFaceIso::calcSubFace
(
    const label& faceI,
    const scalar& cutValue
)
{
    clearStorage();
    const face& f = mesh_.faces()[faceI];
    label inLiquid = 0;
    label firstFullySubmergedPoint = -1;

    // loop over face
    forAll(f, i)
    {
        // pointStatus is f - cutValue
        pointStatus_.append(f_[f[i]] - cutValue);
        if (mag(pointStatus_[i]) < 10 * SMALL)
        {
            pointStatus_[i] = 0;
        }
        if (pointStatus_[i] > 10 * SMALL)
        {
            inLiquid++;
            if (firstFullySubmergedPoint == -1)
            {
                firstFullySubmergedPoint = i;
            }
        }
    }

    if (inLiquid == f.size()) // fluid face
    {
        faceStatus_ = -1;
        subFaceCentre_ = mesh_.faceCentres()[faceI];
        subFaceArea_ = mesh_.faceAreas()[faceI];
        return faceStatus_;
    }
    else if (inLiquid == 0) // gas face
    {
        faceStatus_ = 1;
        subFaceCentre_ = vector::zero;
        subFaceArea_ = vector::zero;
        return faceStatus_;
    }


    cutFace::calcSubFace
    (
        faceI,
        pointStatus_,
        firstFullySubmergedPoint,
        subFacePoints_,
        surfacePoints_,
        faceStatus_,
        subFaceCentre_,
        subFaceArea_
    );

    return faceStatus_;
}

Foam::point Foam::cutFaceIso::subFaceCentre()
{
    return subFaceCentre_;
}

Foam::vector Foam::cutFaceIso::subFaceArea()
{
    return subFaceArea_;
}

Foam::DynamicList<Foam::point>& Foam::cutFaceIso::subFacePoints()
{
    return subFacePoints_;
}

Foam::DynamicList<Foam::point>& Foam::cutFaceIso::surfacePoints()
{
    return surfacePoints_;
}

void Foam::cutFaceIso::clearStorage()
{
    subFaceCentre_ = vector::zero;
    subFaceArea_ = vector::zero;
    subFacePoints_.clear();
    surfacePoints_.clear();
    pointStatus_.clear();
    weight_.clear();
    faceStatus_ = -1;
}

// ************************************************************************* //
