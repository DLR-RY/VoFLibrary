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


#include "gradAlpha.H"
#include "addToRunTimeSelectionTable.H"

#include "fvc.H"
#include "leastSquareGrad.H"



// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reconstruction
{
    defineTypeNameAndDebug(gradAlpha, 0);
    addToRunTimeSelectionTable(reconstructionSchemes,gradAlpha, components);
}
}


void Foam::reconstruction::gradAlpha::gradSurf(const volScalarField& phi)
{
    leastSquareGrad<scalar> lsGrad("polyDegree1",mesh_.geometricD());

    exchangeFields_.setUpCommforZone(interfaceCell_,false);

    Map<vector> mapCC(exchangeFields_.getDatafromOtherProc(interfaceCell_,mesh_.C()));
    Map<scalar> mapPhi(exchangeFields_.getDatafromOtherProc(interfaceCell_,phi));

    DynamicField<vector > cellCentre(100); // should be big enough avoids resizing
    DynamicField<scalar > phiValues(100);

    const labelListList& stencil = exchangeFields_.getStencil();

    forAll(interfaceLabels_, i)
    {
        //if(interfaceCell_[celli])
        //{
        const label celli = interfaceLabels_[i]; 
        cellCentre.clear(); 
        phiValues.clear();

        forAll(stencil[celli],i)
        {
            const label& gblIdx = stencil[celli][i];
            cellCentre.append(exchangeFields_.getValue(mesh_.C(),mapCC,gblIdx));
            phiValues.append(exchangeFields_.getValue(phi,mapPhi,gblIdx));
        }

        cellCentre -= mesh_.C()[celli];
        interfaceNormal_[i] = lsGrad.grad(cellCentre,phiValues);
        // Info << " grad celli " << celli  << " vector " << interfaceNormal_[i] << endl;
        
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reconstruction::gradAlpha::gradAlpha
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U,
    dictionary& dict
)
:
    reconstructionSchemes
    (
        typeName,
        alpha1,
        phi,
        U,
        dict
    ),
    mesh_(alpha1.mesh()),
    interfaceNormal_(fvc::grad(alpha1)),
    vof2IsoTol_(readScalar(modelDict().lookup("vof2IsoTol" ))),
    surfCellTol_(readScalar(modelDict().lookup("surfCellTol" ))),
    exchangeFields_(zoneDistribute::New(mesh_)),
    sIterPLIC_(mesh_,surfCellTol_)
    




{
  writeVTK_ =  readBool( modelDict().lookup("writeVTK" ));
  reconstruct();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reconstruction::gradAlpha::~gradAlpha()
{

}

// * * * * * * * * * * * * * * Protected Access Member Functions  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //
void Foam::reconstruction::gradAlpha::reconstruct()
{
    bool uptodate = alreadyReconstructed();

    if(uptodate)
    {
        return;
    }
    
    if (mesh_.topoChanging())
    {
        // Introduced resizing to cope with changing meshes
        if(interfaceCell_.size() != mesh_.nCells())
        {
            interfaceCell_.resize(mesh_.nCells());
        }
    }
    interfaceCell_ = false;

    interfaceLabels_.clear();
    //interfaceNormal_.clear();
    forAll(alpha1_,celli)
    {
        if(sIterPLIC_.isASurfaceCell(alpha1_[celli]))
        {
            interfaceCell_[celli] = true; // is set to false earlier
            interfaceLabels_.append(celli);
        }
    }
    interfaceNormal_.setSize(interfaceLabels_.size());

    gradSurf(alpha1_);

    forAll(interfaceLabels_, i)
    {
        const label celli = interfaceLabels_[i];
        if(mag(interfaceNormal_[i]) == 0)
        {
            continue;
        } 

        sIterPLIC_.vofCutCell
        (
            celli,
            alpha1_[celli],
            vof2IsoTol_,
            100,
            interfaceNormal_[i]
        );

        if(sIterPLIC_.cellStatus() == 0)
        {

            normal_[celli] = sIterPLIC_.surfaceArea();
            centre_[celli] = sIterPLIC_.surfaceCentre();
            if(mag(normal_[celli]) == 0)
            { 
                normal_[celli] = vector::zero;
                centre_[celli] = vector::zero;
            }

            //interfaceCell_[cellI]=true;
        }
        else
        {
            //interfaceNormal_[i] = vector::zero;
            normal_[celli] = vector::zero;
            centre_[celli] = vector::zero;
        }
    }

}

