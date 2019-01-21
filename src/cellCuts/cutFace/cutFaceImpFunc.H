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
    Foam::cutFaceImpFunc

Description
    Class for cutting a face, faceI, of an fvMesh, mesh_, at its intersection
    with an isosurface defined by the mesh point values f_ and the isovalue,
    cutValue.

    @article{Roenby160405,
        author={Roenby, Johan and Bredmose, Henrik and Jasak, Hrvoje},
        title={A computational method for sharp interface advection},
        volume={3},
        year={2016},
        doi={10.1098/rsos.160405},
        publisher = {The Royal Society},
        journal = {Royal Society Open Science}
    }

Author
    Johan Roenby, DHI, all rights reserved.

SourceFiles
    cutFaceImpFunc.C

\*---------------------------------------------------------------------------*/

#ifndef cutFaceImpFunc_H
#define cutFaceImpFunc_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "OFstream.H"
#include "cutFace.H"
#include "implicitFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class cutFaceImpFuncs Declaration
\*---------------------------------------------------------------------------*/

class cutFaceImpFunc
:
        public cutFace
{
    // Private data

        //- Mesh whose cells and faces to cut at their intersection with an isoface
        const fvMesh& mesh_;

        //- Isofunction values at mesh points. f_size() = mesh_.nPoints().
        scalarField& f_;

        // implicit function
        const implicitFunctions& func_;

        //- Storage for centre of subface
        point subFaceCentre_;

        //- Storage for area vector of subface
        vector subFaceArea_;

        //- Storage for subFacePoints
        DynamicList<point> subFacePoints_;

        //- Storage for subFacePoints
        DynamicList<point> surfacePoints_;

        //- Storage for pointStatus_ cuts the cell at 0
        DynamicList<scalar> pointStatus_;

        //- Storage of the edge weight
        DynamicList<scalar> weight_;


        //- A face status label taking one of the values:
        //
        //  - -1: face is fully below the isosurface
        //  -  0: face is cut, i.e. has values larger and smaller than isoValue_
        //  - +1: face is fully above the isosurface
        label faceStatus_;


public:

    // Constructors

        //- Construct from fvMesh and a scalarField
        //  Length of scalarField should equal number of mesh points
        cutFaceImpFunc(const fvMesh& mesh,scalarField& f,const implicitFunctions& func);


    // Member functions
        scalar bisection(point p0, point p1);

        //- Calculate cut points along edges of faceI
        label calcSubFace
        (
            const label& faceI,
            const scalar& cutValue
        );

        //- returns centre of cutted face
        point subFaceCentre();

        //- returns area vector of cutted face
        vector subFaceArea();

        //- returns the cut edge of the cutted face
        DynamicList<point>& subFacePoints();

        //- returns point of the face in sorted of cutted face
        DynamicList<point>& surfacePoints();

        //- resets internal variables
        void clearStorage();
        

        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
