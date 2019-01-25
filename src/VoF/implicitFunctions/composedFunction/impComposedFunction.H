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
    Foam::impComposedFunction

Description
    Handles multiple implict functions and offers multiple ways to combine
    them

Usage
    Example of function object partial specification:
     \verbatim
    function composedFunction;
    mode minDist;
    // following mode are avaible:
    // "add" "subtract" "minDist" "intersect"
    composedFunctions
    {
        plane
        {
            function plane;
            centre (0 1. 0);
            normal (0 1 0);
            
        }

        sphere
        {
            function sphere;
            radius 0.4;
            centre (0.5 1.5 0.5);
            scale 1;
        }

        sphere2
        {
            function sphere;
            radius 0.4;
            centre (0.5 0.5 0.5);
            scale -1;
        }
    }
    \endverbatim




Author
    Henning Scheufler, DLR, all rights reserved.

SourceFiles
    impComposedFunction.C

\*---------------------------------------------------------------------------*/

#ifndef impComposedFunction_H
#define impComposedFunction_H

#include "implicitFunctions.H"
#include "mathematicalConstants.H"
#include "cartesianCS.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace implicitFunction
{
// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class impComposedFunction Declaration
\*---------------------------------------------------------------------------*/

class impComposedFunction
:
    public implicitFunctions
{
private:

    // Private Member Data
             //- Enumeration defining the valid actions
    enum setMode
    {
        ADD,                
        SUBTRACT,           
        MINDIST,
        INTERSECT 
    };

    //- The setActions text representations
    static const Enum<setMode> modeNames;

    // stores the functions
    PtrList<implicitFunctions> functions_;

    setMode mode_;

    // needs for finding the closest function. 
    // avoid creation every call
    scalarField values_;

    label selectFunction(const scalarField& values_);


        

    // Private Member Functions

        //- Disallow default bitwise copy construct
        impComposedFunction(const impComposedFunction&);

        //- Disallow default bitwise assignment
        void operator=(const impComposedFunction&);


public:

    //- Runtime type information
    TypeName("composedFunction");


    // Constructors

        // //- Construct from components
        // impComposedFunction
        // (
        //     const vector centre,
        //     const vector normal,
        //     const vector tangent,
        //     const scalar k1,
        //     const scalar k2
        // );

        //- Construct from dictionary (used by implicitFunctions)
        impComposedFunction
        (
            const dictionary& dict
        );

    //- Destructor
    virtual ~impComposedFunction();

    virtual scalar value(const vector p);

    virtual vector grad(const vector p);

    virtual scalar distanceToSurfaces(const vector p);


    // Member Functions


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace implicitFunction
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
