/*---------------------------------------------------------------------------*\
    Modified work | Copyright (c) 2017-2019, German Aerospace Center (DLR)
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
    Foam::isoAlpha

Description
    Reconstructs an interface (centre and normal vectors) consisting of isosurfaces
    to match the internal fluid distribution in cells. The point values (ap_)
    are estimated by inverse distance interpolation of the VoF (alpha) field.

    Reference:
        \verbatim
            Roenby, J., Bredmose, H. and Jasak, H. (2016).
            A computational method for sharp interface advection
            Royal Society Open Science, 3
            doi 10.1098/rsos.160405
        \endverbatim

    Original code supplied by Johan Roenby, DHI (2016)

Author
    Johan Roenby, DHI, all rights reserved.
    Modified Henning Scheufler, DLR

SourceFiles
    isoAlpha.C

\*---------------------------------------------------------------------------*/

#ifndef isoAlpha_H
#define isoAlpha_H

#include "autoPtr.H"
#include "dimensionedScalar.H"
#include "reconstructionSchemes.H"
#include "typeInfo.H"
#include "volFields.H"

#include "surfaceIteratorIso.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reconstruction
{

/*---------------------------------------------------------------------------*\
                Class isoAlpha Declaration
\*---------------------------------------------------------------------------*/

class isoAlpha
:
    public reconstructionSchemes
{

  private:
        //- Reference to mesh
        const fvMesh& mesh_;

        //- Interpolation object from cell centres to points
        volPointInterpolation vpi_;

        //- VOF field interpolated to mesh points
        scalarField ap_;

        // Switches and tolerances. Tolerances need to go into toleranceSwitches

        //- Tolerance for search of isoFace giving specified VOF value
        scalar vof2IsoTol_;

        //- Tolerance for marking of surface cells:
        //  Those with surfCellTol_ < alpha1 < 1 - surfCellTol_
        scalar surfCellTol_;

        //- surfaceIterator finds the isovalue for specified VOF value
        surfaceIteratorIso sIterIso_;


        // Private Member Functions

        //- Disallow copy construct
        isoAlpha(const isoAlpha&);

        //- Disallow default bitwise assignment
        void operator=(const isoAlpha&);

  public:
        //- Runtime type information
        TypeName("isoAlpha");

        // Constructors

        //- Construct from components
        isoAlpha
        (
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U,
            dictionary& dict
        );

        //- Destructor
        virtual ~isoAlpha();

        //- reconstructs the interface
        virtual void reconstruct();

        virtual void mapAlphaField()
        {
            // do nothing;
        }
        
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End reconstructionSchemes Foam
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
