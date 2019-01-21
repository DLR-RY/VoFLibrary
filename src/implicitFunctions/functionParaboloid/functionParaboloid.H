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
    Foam::functionParaboloid

Description

Author
    Henning Scheufler, DLR

SourceFiles
    functionParaboloid.C

\*---------------------------------------------------------------------------*/

#ifndef functionParaboloid_H
#define functionParaboloid_H

#include "implicitFunctions.H"
#include "mathematicalConstants.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class functionParaboloid Declaration
\*---------------------------------------------------------------------------*/

class functionParaboloid
:
    public implicitFunctions
{
private:

    // Private Member Data

        //- ax^2 + bx*y + cy^2
        vector coeffs_;



    // Private Member Functions

        //- Disallow default bitwise copy construct
        functionParaboloid(const functionParaboloid&);

        //- Disallow default bitwise assignment
        void operator=(const functionParaboloid&);


public:

    //- Runtime type information
    TypeName("functionParaboloid");


    // Constructors

        //- Construct from components
        functionParaboloid(const vector coeffs);

        //- Construct from dictionary (used by implicitFunctions)
        functionParaboloid
        (
            const dictionary& dict
        );

    //- Destructor
    virtual ~functionParaboloid();

    virtual scalar value(const vector p) const
    {
          return coeffs_.x()*sqr(p.x()) + coeffs_.y()*p.x()*p.y()  + coeffs_.z()*sqr(p.y()) - p.z();
    }


    // Member Functions


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
