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
    Foam::implicitFunctions

Description
    base class for implicit functions

Author
    Henning Scheufler, DLR, all rights reserved.

SourceFiles
    implicitFunctions.C

\*---------------------------------------------------------------------------*/

#ifndef implicitFunctions_H
#define implicitFunctions_H


#include "runTimeSelectionTables.H"
#include "objectRegistry.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class implicitFunctions Declaration
\*---------------------------------------------------------------------------*/

class implicitFunctions
{
    // Private data



    // Private Member Functions

        //- Disallow default bitwise copy construct
        implicitFunctions(const implicitFunctions&);

        //- Disallow default bitwise assignment
        void operator=(const implicitFunctions&);


public:

    //- Runtime type information
    TypeName("implicitFunctions");

    // Declare run-time constructor selection table

        // For the dictionary constructor
        declareRunTimeSelectionTable
        (
            autoPtr,
            implicitFunctions,
            dict,
            (
                const dictionary& dict
            ),
            (dict)
        );


 


    // Constructors

        implicitFunctions();

     // Selectors

        //- Return a reference to the selected implicitFunctions
        static autoPtr<implicitFunctions> New
        (
            const word& implicitFunctionsType,
            const dictionary& dict
        );


    //- Destructor
    virtual ~implicitFunctions();


    // Member Functions

        //- Is object global
        virtual scalar value(const vector p) 
        {
            return GREAT;
        }

        virtual vector grad(const vector p) 
        {
            return vector(GREAT,GREAT,GREAT);
        }

        virtual scalar distanceToSurfaces(const vector p) 
        {
            return GREAT;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
