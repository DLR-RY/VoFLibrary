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

#include "cutFacePLIC.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cutFacePLIC::cutFacePLIC(const fvMesh& mesh)
    : cutFace(mesh),
      mesh_(mesh),
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

Foam::label Foam::cutFacePLIC::calcSubFace
(
    const label& faceI,
    const vector& normal,
    const vector& base
 )
{
    clearStorage();

    const face& f = mesh_.faces()[faceI];
    label inLiquid = 0;
    label firstFullySubmergedPoint = -1;

    // loop face
    forAll(f, i)
    {
        // pointStatus is the distance to the plane
        scalar value = (mesh_.points()[f[i]] - base) & normal;
        if (mag(value) < SMALL)
        {
            value = 0;
        }

        pointStatus_.append(value);
        if (pointStatus_[i] > 0)
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


Foam::point Foam::cutFacePLIC::subFaceCentre()
{
    return subFaceCentre_;
}

Foam::vector Foam::cutFacePLIC::subFaceArea()
{
    return subFaceArea_;
}

Foam::DynamicList<Foam::point>& Foam::cutFacePLIC::subFacePoints()
{
    return subFacePoints_;
}

Foam::DynamicList<Foam::point>& Foam::cutFacePLIC::surfacePoints()
{
    return surfacePoints_;
}

void Foam::cutFacePLIC::clearStorage()
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
