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
    Foam::implicitFunction::impEllipsoid

Description
    creates an ellipsoid

Author
    Henning Scheufler, DLR, all rights reserved.

SourceFiles
    impEllipsoid.C

\*---------------------------------------------------------------------------*/

#ifndef impEllipsoid_H
#define impEllipsoid_H

#include "implicitFunctions.H"
#include "mathematicalConstants.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace implicitFunction
{
// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class impEllipsoid Declaration
\*---------------------------------------------------------------------------*/

class impEllipsoid
:
    public implicitFunctions
{
private:

    // Private Member Data

        //- Centre point
        vector semiAxis_;

        //- Centre point
        vector centre_;



    // Private Member Functions

        //- Disallow default bitwise copy construct
        impEllipsoid(const impEllipsoid&);

        //- Disallow default bitwise assignment
        void operator=(const impEllipsoid&);


public:

    //- Runtime type information
    TypeName("ellipsoid");


    // Constructors

        //- Construct from components
        impEllipsoid(const vector semiAxis);

        //- Construct from dictionary (used by implicitFunctions)
        impEllipsoid
        (
            const dictionary& dict
        );

    //- Destructor
    virtual ~impEllipsoid();

    virtual scalar value(const vector p) 
    {
          return -sqrt(sqr((p.x()-centre_.x())/semiAxis_.x())+sqr((p.y()-centre_.y())/semiAxis_.y())+sqr((p.z()-centre_.z())/semiAxis_.z())) + 1;
    }

    virtual vector grad(const vector p) 
    {
        // normal_ has the length of one
        return vector
        (
            2*(p.x()-centre_.x())/sqr(semiAxis_.x()),
            2*(p.y()-centre_.y())/sqr(semiAxis_.y()),
            2*(p.z()-centre_.z())/sqr(semiAxis_.z())
        );
    }

    virtual scalar distanceToSurfaces(const vector p) 
    {
        // not yet implemented
        return GREAT;
    }


    // Member Functions


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace implicitFunction
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
