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
    Foam::functionSin

Description

Author
    Henning Scheufler, DLR

SourceFiles
    functionSin.C

\*---------------------------------------------------------------------------*/

#ifndef functionSin_H
#define functionSin_H

#include "implicitFunctions.H"
#include "mathematicalConstants.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class functionSin Declaration
\*---------------------------------------------------------------------------*/

class functionSin
:
    public implicitFunctions
{
private:

    // Private Member Data

        //- Centre point
        const scalar period_;

        //- Radius
        const scalar phase_;

        const scalar amplitude_;

        vector up_;

        vector direction_;

        vector centre_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        functionSin(const functionSin&);

        //- Disallow default bitwise assignment
        void operator=(const functionSin&);


public:

    //- Runtime type information
    TypeName("functionSin");


    // Constructors

        //- Construct from components
        functionSin
        (
            const scalar period,
            const scalar phase,
            const scalar amplitude,
            vector direction,
            vector up,
            vector centre
        );

        //- Construct from dictionary (used by implicitFunctions)
        functionSin
        (
            const dictionary& dict
        );

    //- Destructor
    virtual ~functionSin();

    virtual scalar value(const vector p) const
    {

        scalar x = (p - centre_) & direction_;
        scalar z = (p - centre_) & up_;

        return z - amplitude_*Foam::sin(2*constant::mathematical::pi*x/period_ + phase_);
    }


    // Member Functions


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
