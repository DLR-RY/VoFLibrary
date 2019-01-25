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

\*---------------------------------------------------------------------------*/

#include "impSin.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace implicitFunction
    {
        defineTypeNameAndDebug(impSin, 0);
        addToRunTimeSelectionTable(implicitFunctions, impSin, dict);
    }

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::implicitFunction::impSin::impSin
(
    const scalar period,
    const scalar phase,
    const scalar amplitude,
    const vector direction,
    const vector up,
    const vector centre
)
:
    period_(period),
    phase_(phase),
    amplitude_(amplitude),
    up_(up),
    direction_(direction),
    centre_(centre)
{

}


Foam::implicitFunction::impSin::impSin
(
    const dictionary& dict
)
:
    period_(readScalar(dict.lookup("period"))),
    phase_(readScalar(dict.lookup("phase"))),
    amplitude_(readScalar(dict.lookup("amplitude"))),
    up_(dict.lookup("up")),
    direction_(dict.lookup("direction")),
    centre_(dict.lookup("centre"))
{
    direction_ /= mag(direction_);
    up_ /= mag(up_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::implicitFunction::impSin::~impSin()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
