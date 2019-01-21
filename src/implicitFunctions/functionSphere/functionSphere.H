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
    Foam::functionSphere

Description

Author
    Henning Scheufler, DLR

SourceFiles
    functionSphere.C

\*---------------------------------------------------------------------------*/

#ifndef functionSphere_H
#define functionSphere_H

#include "implicitFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class functionSphere Declaration
\*---------------------------------------------------------------------------*/

class functionSphere
:
    public implicitFunctions
{
private:

    // Private Member Data

        //- Centre point
        const point centre_;

        //- Radius
        const scalar radius_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        functionSphere(const functionSphere&);

        //- Disallow default bitwise assignment
        void operator=(const functionSphere&);


public:

    //- Runtime type information
    TypeName("functionSphere");


    // Constructors

        //- Construct from components
        functionSphere( const point&, const scalar radius);

        //- Construct from dictionary (used by implicitFunctions)
        functionSphere
        (
            const dictionary& dict
        );

    //- Destructor
    virtual ~functionSphere();

    // Member Functions

    virtual scalar value(const vector p) const
    {
          return -mag(p-centre_)+radius_;
    }

    virtual vector grad(const vector p) const
    {
        return centre_-p;
    }

    virtual scalar distanceToSurfaces(const vector p) const
    {
        return mag(mag(p-centre_)-radius_);
    }





};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
