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
    Foam::gradAlpha

Description
    Reconstructs an interface (centre and normal vector) consisting of planes
    to match the internal fluid distribution in cells. The interface normals
    are estimated by least square gradient scheme on the VoF Field (alpha).

Author
    Henning Scheufler, DLR
 
SourceFiles
    gradAlpha.C

\*---------------------------------------------------------------------------*/

#ifndef gradAlpha_H
#define gradAlpha_H

#include "typeInfo.H"
#include "reconstructionSchemes.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"


#include "surfaceIteratorPLIC.H"
#include "heightFunc.H"




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reconstruction
{

/*---------------------------------------------------------------------------*\
                Class gradAlpha Declaration
\*---------------------------------------------------------------------------*/

class gradAlpha
    : public reconstructionSchemes
{

private:

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Interpolation object from cell centres to points
        volVectorField interfaceNormal_;




    // Switches and tolerances. Tolerances need to go into toleranceSwitches


        //- Tolerance for search of isoFace giving specified VOF value
        scalar vof2IsoTol_;

        //- Tolerance for marking of surface cells:
        //  Those with surfCellTol_ < alpha1 < 1 - surfCellTol_
        scalar surfCellTol_;

        // calculates the RDF function? why do I need this

        //- surfaceIterator finds the plane centre for specified VOF value
        surfaceIteratorPLIC sIterPLIC_;

        // calculates the RDF function
        heightFunc heightFunc_;


private:

    // Private Member Functions

        //- Disallow copy construct
        gradAlpha(const gradAlpha&);

        //- Disallow default bitwise assignment
        void operator=(const gradAlpha&);



public:

    //- Runtime type information
    TypeName("gradAlpha");


    // Constructors

        //- Construct from components
        gradAlpha
        (
                const volScalarField& alpha1,
                const surfaceScalarField& phi,
                const volVectorField& U,
                dictionary& dict
        );


    //- Destructor
    virtual ~gradAlpha();

    // reconstruct  interface
    virtual void reconstruct();


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End reconstructionSchemes Foam
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
