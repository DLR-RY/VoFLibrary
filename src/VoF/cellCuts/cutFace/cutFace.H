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
    Foam::cutFace

Description
    Base class for cutting a face, faceI, of an fvMesh, mesh_, at its
    intersections

Author
    Henning Scheufler, DLR, all rights reserved.
    Johan Roenby, DHI, all rights reserved.


SourceFiles
    cutFace.C

\*---------------------------------------------------------------------------*/

#ifndef cutFace_H
#define cutFace_H

#include "OFstream.H"
#include "fvMesh.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class cutFaces Declaration
\*---------------------------------------------------------------------------*/

class cutFace
{
  private:
    // Private data

    //- Reference to mesh
    const fvMesh& mesh_;

  protected:

    //- Calculate cut points along edges of face with the
    // pointStatus and the face label
    //  Returns the face status, where:

    //  -1: faceStatus is fully below the isosurface
    //   0: faceStatus is cut, i.e. has values larger and smaller than isoValue
    //  +1: faceStatus is fully above the isosurface
    void calcSubFace
    (
        const label& faceI,
        const scalarList& pointStatus,
        label firstFullySubmergedPoint,
        DynamicList<point>& subFacePoints,
        DynamicList<point>& surfacePoints,
        label& faceStatus,
        vector& subFaceCentre,
        vector& subFaceArea
     );

    //  -1: faceStatus is fully below the isosurface
    //   0: faceStatus is cut, i.e. has values larger and smaller than isoValue
    //  +1: faceStatus is fully above the isosurface
    void calcSubFace
    (
        const label& faceI,
        const scalarList& pointStatus,
        const scalarList& weights,
        label firstFullySubmergedPoint,
        DynamicList<point>& subFacePoints,
        DynamicList<point>& surfacePoints,
        label& faceStatus,
        vector& subFaceCentre,
        vector& subFaceArea
     );

    // calculates centre and normal of the face
    void calcSubFaceCentreAndArea
    (
        DynamicList<point>& subFacePoints,
        vector& subFaceCentre,
        vector& subFaceArea
    );

    //- Calculate cut points along edges of face with
    // pointStatus, pointfield and face
    //-  Returns the face status, where:

    // -  -1: faceStatus is fully below the isosurface
    // -   0: faceStatus is cut, i.e. has values larger and smaller than isoValue
    // -  +1: faceStatus is fully above the isosurface
    void calcSubFace
    (
        const face& f,
        const pointField& points,
        const scalarList& pointStatus,
        label firstFullySubmergedPoint,
        DynamicList<point>& subFacePoints,
        DynamicList<point>& surfacePoints,
        label& faceStatus,
        vector& subFaceCentre,
        vector& subFaceArea
    );

  public:
    // Constructors

    //- Construct from fvMesh
    cutFace(const fvMesh& mesh);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
