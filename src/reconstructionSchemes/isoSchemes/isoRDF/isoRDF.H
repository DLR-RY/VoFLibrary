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
    Foam::isoRDF

Description
    Reconstructs an interface (centre and normal vectors) consisting of isosurfaces
    to match the internal fluid distribution in cells. The point values (ap_)
    are estimated by least square interpolation of the RDF function (height).

Author
    Henning Scheufler, DLR
 
SourceFiles
    isoRDF.C

\*---------------------------------------------------------------------------*/

#ifndef isoRDF_H
#define isoRDF_H

#include "autoPtr.H"
#include "dimensionedScalar.H"
#include "reconstructionSchemes.H"
#include "typeInfo.H"
#include "volFields.H"

#include "surfaceIteratorIso.H"
#include "volPointInterpolation.H"

#include "heightFunc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reconstruction
{

/*---------------------------------------------------------------------------*\
                Class isoRDF Declaration
\*---------------------------------------------------------------------------*/

class isoRDF : public reconstructionSchemes
{

  private:
    //- Reference to mesh
    const fvMesh& mesh_;

    //- Interpolation object from cell centres to points
    volPointInterpolation vpi_;

    //- VOF field interpolated to mesh points
    scalarField ap_;

    //- Interpolation object from cell centres to points
    volScalarField height_;

    // Switches and tolerances. Tolerances need to go into toleranceSwitches

    //- Tolerance for search of isoFace giving specified VOF value
    scalar vof2IsoTol_;

    //- Tolerance for marking of surface cells:
    //  Those with surfCellTol_ < alpha1 < 1 - surfCellTol_
    scalar surfCellTol_;

    //  tolerance
    scalar tol_;

    // number of times that the interface ist reconstructed has to be bigger than 2
    label iteration_;

    // calculates the RDF function
    heightFunc heightFunc_;

    //- surfaceIterator finds the isovalue for specified VOF value
    surfaceIteratorIso sIterIso_;

  private:
    // Private Member Functions

    //- Disallow copy construct
    isoRDF(const isoRDF&);

    //- Disallow default bitwise assignment
    void operator=(const isoRDF&);

  public:
    //- Runtime type information
    TypeName("isoRDF");

    // Constructors

    //- Construct from components
    isoRDF
    (
        const volScalarField& alpha1,
        const surfaceScalarField& phi,
        const volVectorField& U,
        dictionary& dict
    );

    //- Destructor
    virtual ~isoRDF();

    //- constructs the interface
    virtual void reconstruct();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End reconstructionSchemes Foam
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
