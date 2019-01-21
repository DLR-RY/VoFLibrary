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
    Foam::functionDisc

Description

Author
    Henning Scheufler, DLR

SourceFiles
    functionDisc.C

\*---------------------------------------------------------------------------*/

#ifndef functionDisc_H
#define functionDisc_H

#include "implicitFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class functionDisc Declaration
\*---------------------------------------------------------------------------*/

class functionDisc
:
    public implicitFunctions
{
private:

    // Private Member Data

        //- Centre point
        const point centre_;

        //- Radius
        const scalar radius_;

        vector direction_;

        tensor project_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        functionDisc(const functionDisc&);

        //- Disallow default bitwise assignment
        void operator=(const functionDisc&);


public:

    //- Runtime type information
    TypeName("functionDisc");


    // Constructors

        //- Construct from components
        functionDisc( const point&, const scalar radius,const vector direction);

        //- Construct from dictionary (used by implicitFunctions)
        functionDisc
        (
            const dictionary& dict
        );

    //- Destructor
    virtual ~functionDisc();

    virtual scalar value(const vector p) const
    {
          return -mag(project_ & (p-centre_))+radius_;
    }

    virtual vector grad(const vector p) const
    {
        return -(project_ & (p-centre_));
    }

    virtual scalar distanceToSurfaces(const vector p) const
    {
        return mag(mag(project_ & (p-centre_))-radius_);
    }


    // Member Functions


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
