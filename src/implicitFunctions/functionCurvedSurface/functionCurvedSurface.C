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

#include "functionCurvedSurface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(functionCurvedSurface, 0);
addToRunTimeSelectionTable(implicitFunctions, functionCurvedSurface, dict);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionCurvedSurface::functionCurvedSurface
(
    const vector centre,
    const vector normal,
    const vector tangent,
    const scalar k1,
    const scalar k2
)
:
  centre_(centre),
  normal_(normal),
  tangent_(tangent),
  k1_(k1),
  k2_(k2)
{
    coordSys_ = cartesianCS(word("coordSys"),centre_,normal_,tangent_);
}


Foam::functionCurvedSurface::functionCurvedSurface
(
    const dictionary& dict
)
:
    centre_(dict.lookup("centre")),
    normal_(dict.lookup("normal")),
    tangent_(dict.lookup("tangent")),
    k1_(readScalar(dict.lookup("k1"))),
    k2_(readScalar(dict.lookup("k2")))
{
    coordSys_ = cartesianCS(word("coordSys"),centre_,normal_,tangent_);

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionCurvedSurface::~functionCurvedSurface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
