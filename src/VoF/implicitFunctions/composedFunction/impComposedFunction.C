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

#include "impComposedFunction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace implicitFunction
    {
        defineTypeNameAndDebug(impComposedFunction, 0);
        addToRunTimeSelectionTable(implicitFunctions, impComposedFunction, dict);
    }

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

 const Foam::Enum
 <
    Foam::implicitFunction::impComposedFunction::setMode
 >
 Foam::implicitFunction::impComposedFunction::modeNames
 ({
     { setMode::ADD, "add" },
     { setMode::SUBTRACT, "subtract" },
     { setMode::MINDIST, "minDist" },
     { setMode::INTERSECT, "intersect" }, //,
 });

Foam::label Foam::implicitFunction::impComposedFunction::selectFunction(const scalarField& values)
{
    switch (mode_)
    {
        case setMode::MINDIST:
        {
            scalarField absVal = mag(values);
            return findMin(absVal);
        }
        case setMode::ADD:
        {
            return findMax(values);
        }
        case setMode::SUBTRACT:
        {
            label idx = findMin(values,1); // start at the second entry
            if(values[idx] < values[0] && pos(values[0]))
            {
                return idx;
            }
            else
            {
                return 0;
            }
        }
        case setMode::INTERSECT:
        {
            return findMin(values);
        }
        default:
        {
            FatalErrorInFunction
                << "This mode is not supported  only " << nl
                << "Supported modes are: " << nl
                << "minDist" << nl
                << "add" << nl
                << "subtract" << nl
                << abort(FatalError);
            return -1;
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //



Foam::implicitFunction::impComposedFunction::impComposedFunction
(
    const dictionary& dict
)
:
    functions_(dict.subDict("composedFunctions").size()),
    //mode_(dict.subDict("impComposedFunctions").get<word>("mode")),
    mode_(modeNames.get("mode",dict)),
    values_(dict.subDict("composedFunctions").size())
{
    dictionary funcDict = dict.subDict("composedFunctions");
    label funcI = 0;

    forAllConstIter(dictionary,funcDict,iter)
    {
        const word& key = iter().keyword();

        if (!funcDict.isDict(key))
        {
            FatalErrorInFunction
                << "Found non-dictionary entry " << iter()
                << " in top-level dictionary " << funcDict
                << exit(FatalError);
        }

        const dictionary& compFuncDict = funcDict.subDict(key);

        functions_.set
        (
            funcI,
            implicitFunctions::New
            (
                word(compFuncDict.lookup("function")),
                compFuncDict
            )
        );

        funcI++;

    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::implicitFunction::impComposedFunction::~impComposedFunction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::scalar Foam::implicitFunction::impComposedFunction::value(const vector p)
{

    forAll(values_,i)
    {
        // values_[i] = mag(functions_[i].value(p));
        values_[i] = (functions_[i].value(p));
    }
    label idx = selectFunction(values_);

    // return functions_[minIdx].value(p);
    return values_[idx];
}

Foam::vector Foam::implicitFunction::impComposedFunction::grad(const vector p)
{

    forAll(values_,i)
    {
        values_[i] = mag(functions_[i].value(p));
    }
    label minIdx = findMin(values_);

    return functions_[minIdx].grad(p);
}

Foam::scalar Foam::implicitFunction::impComposedFunction::distanceToSurfaces(const vector p)
{

    forAll(values_,i)
    {
        values_[i] = mag(functions_[i].value(p));
    }
    label minIdx = findMin(values_);

    return functions_[minIdx].distanceToSurfaces(p);
}

// ************************************************************************* //
