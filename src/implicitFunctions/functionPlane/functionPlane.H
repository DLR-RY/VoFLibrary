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
    Foam::functionPlane

Description

Author
    Henning Scheufler, DLR

SourceFiles
    functionPlane.C

\*---------------------------------------------------------------------------*/

#ifndef functionPlane_H
#define functionPlane_H

#include "implicitFunctions.H"
#include "mathematicalConstants.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class functionPlane Declaration
\*---------------------------------------------------------------------------*/

class functionPlane
:
    public implicitFunctions
{
private:

    // Private Member Data

        //- ax^2 + bx*y + cy^2
    // reference point of the plane
        vector centre_;
    // normal vector of the plane
        vector normal_;




    // Private Member Functions

        //- Disallow default bitwise copy construct
        functionPlane(const functionPlane&);

        //- Disallow default bitwise assignment
        void operator=(const functionPlane&);


public:

    //- Runtime type information
    TypeName("functionPlane");


    // Constructors

        //- Construct from components
        functionPlane(const vector centre,const vector normal);

        //- Construct from dictionary (used by implicitFunctions)
        functionPlane
        (
            const dictionary& dict
        );

    //- Destructor
    virtual ~functionPlane();

    virtual scalar value(const vector p) const
    {
        return normal_ & (centre_ - p);
    }

    virtual vector grad(const vector p) const
    {
        // normal_ has the length of one
        return -normal_;
    }

    virtual scalar distanceToSurfaces(const vector p) const
    {
        // normal_ has the length of one
        return mag((p - centre_) & normal_);
    }


    // Member Functions


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
