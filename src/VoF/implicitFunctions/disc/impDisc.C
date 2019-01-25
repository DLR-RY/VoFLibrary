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

#include "impDisc.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace implicitFunction
    {
        defineTypeNameAndDebug(impDisc, 0);
        addToRunTimeSelectionTable(implicitFunctions, impDisc, dict);
    }

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::implicitFunction::impDisc::impDisc
(
    const point& centre,
    const scalar radius,
    const scalar scale,
    const vector direction
)
:
    centre_(centre),
    radius_(radius),
    scale_(scale),
    direction_(direction)
{
   direction_ /= mag(direction_);
   project_ = tensor::I - direction_ * direction_; // outer product
}


Foam::implicitFunction::impDisc::impDisc
(
    const dictionary& dict
)
:
    centre_(dict.lookup("centre")),
    radius_(readScalar(dict.lookup("radius"))),
    scale_(dict.lookupOrDefault<scalar>("scale",1)),
    direction_(dict.lookup("direction"))
{
    direction_ /= mag(direction_);
    project_ = tensor::I - (direction_ * direction_); // outer product
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::implicitFunction::impDisc::~impDisc()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
