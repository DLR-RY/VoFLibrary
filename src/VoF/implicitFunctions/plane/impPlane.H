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
    Foam::implicitFunction::impPlane

Description
    creates a plane

Author
    Henning Scheufler, DLR, all rights reserved.

SourceFiles
    impPlane.C

\*---------------------------------------------------------------------------*/

#ifndef impPlane_H
#define impPlane_H

#include "implicitFunctions.H"
#include "mathematicalConstants.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace implicitFunction
{
// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class impPlane Declaration
\*---------------------------------------------------------------------------*/

class impPlane
:
    public implicitFunctions
{
private:

    // Private Member Data

        //- ax^2 + bx*y + cy^2
    // reference point of the impPlane
        vector centre_;
    // normal vector of the impPlane
        vector normal_;




    // Private Member Functions

        //- Disallow default bitwise copy construct
        impPlane(const impPlane&);

        //- Disallow default bitwise assignment
        void operator=(const impPlane&);


public:

    //- Runtime type information
    TypeName("plane");


    // Constructors

        //- Construct from components
        impPlane(const vector centre,const vector normal);

        //- Construct from dictionary (used by implicitFunctions)
        impPlane
        (
            const dictionary& dict
        );

    //- Destructor
    virtual ~impPlane();

    virtual scalar value(const vector p) 
    {
        return normal_ & (centre_ - p);
    }

    virtual vector grad(const vector p) 
    {
        // normal_ has the length of one
        return -normal_;
    }

    virtual scalar distanceToSurfaces(const vector p) 
    {
        // normal_ has the length of one
        return mag((p - centre_) & normal_);
    }


    // Member Functions


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace implicitFunction
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
