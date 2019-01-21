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
    Foam::functionEllipsoid

Description

Author
    Henning Scheufler, DLR

SourceFiles
    functionEllipsoid.C

\*---------------------------------------------------------------------------*/

#ifndef functionEllipsoid_H
#define functionEllipsoid_H

#include "implicitFunctions.H"
#include "mathematicalConstants.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class functionEllipsoid Declaration
\*---------------------------------------------------------------------------*/

class functionEllipsoid
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
        functionEllipsoid(const functionEllipsoid&);

        //- Disallow default bitwise assignment
        void operator=(const functionEllipsoid&);


public:

    //- Runtime type information
    TypeName("functionEllipsoid");


    // Constructors

        //- Construct from components
        functionEllipsoid(const vector semiAxis);

        //- Construct from dictionary (used by implicitFunctions)
        functionEllipsoid
        (
            const dictionary& dict
        );

    //- Destructor
    virtual ~functionEllipsoid();

    virtual scalar value(const vector p) const
    {
          return -sqrt(sqr((p.x()-centre_.x())/semiAxis_.x())+sqr((p.y()-centre_.y())/semiAxis_.y())+sqr((p.z()-centre_.z())/semiAxis_.z())) + 1;
    }

    virtual vector grad(const vector p) const
    {
        // normal_ has the length of one
        return vector
        (
            2*(p.x()-centre_.x())/sqr(semiAxis_.x()),
            2*(p.y()-centre_.y())/sqr(semiAxis_.y()),
            2*(p.z()-centre_.z())/sqr(semiAxis_.z())
        );
    }

    virtual scalar distanceToSurfaces(const vector p) const
    {
        // not yet implemented
        return GREAT;
    }


    // Member Functions


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
