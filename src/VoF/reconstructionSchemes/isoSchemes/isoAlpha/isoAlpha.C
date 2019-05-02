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
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>..

\*---------------------------------------------------------------------------*/


#include "isoAlpha.H"
#include "addToRunTimeSelectionTable.H"
#include "cutCellPLIC.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reconstruction
{
    defineTypeNameAndDebug(isoAlpha, 0);
    addToRunTimeSelectionTable(reconstructionSchemes,isoAlpha, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reconstruction::isoAlpha::isoAlpha
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
    // Interpolation data

    vpi_(mesh_),
    ap_(mesh_.nPoints()),


    // Tolerances and solution controls
//    vof2IsoTol_(modelDict().lookupOrAddDefault<scalar>("vof2IsoTol", 1e-8,false,false)),
//    surfCellTol_(modelDict().lookupOrAddDefault<scalar>("surfCellTol", 1e-8,false,false))
    vof2IsoTol_(readScalar(modelDict().lookup("vof2IsoTol" ))),
    surfCellTol_(readScalar(modelDict().lookup("surfCellTol" ))),
    sIterIso_(mesh_,ap_,surfCellTol_)




{
  writeVTK_ =  readBool( modelDict().lookup("writeVTK" ));
  reconstruct();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reconstruction::isoAlpha::~isoAlpha()
{}

// * * * * * * * * * * * * * * Protected Access Member Functions  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //
void Foam::reconstruction::isoAlpha::reconstruct()
{
    bool uptodate = alreadyReconstructed();

    if(uptodate)
    {
        return;
    }
    
    // Interpolating alpha1 cell centre values to mesh points (vertices)
    if (mesh_.topoChanging())
    {
        // Introduced resizing to cope with changing meshes
        if(ap_.size() != mesh_.nPoints())
        {
            ap_.resize(mesh_.nPoints());

        }
        if(interfaceCell_.size() != mesh_.nCells())
        {
            interfaceCell_.resize(mesh_.nCells());
        }
    }
    ap_ = vpi_.interpolate(alpha1_);

    DynamicList< List<point> > facePts;



//    vector avgCentre= vector::zero;
    interfaceLabels_.clear();



    forAll(alpha1_,cellI)
    {
    //
        if(sIterIso_.isASurfaceCell(alpha1_[cellI]))
        {
            interfaceLabels_.append(cellI);


            sIterIso_.vofCutCell
            (
                    cellI,
                    alpha1_[cellI],
                    vof2IsoTol_,
                    100
            );

            if(sIterIso_.cellStatus() == 0)
            {

                normal_[cellI] = sIterIso_.surfaceArea();
                centre_[cellI] = sIterIso_.surfaceCentre();
                if(mag(normal_[cellI]) != 0)
                {
                  interfaceCell_[cellI]=true;
                //  facePts.append(sIterIso_.facePoints());
                  if (writeVTK_ & mesh_.time().writeTime())
                  {
                          facePts.append(sIterIso_.facePoints());
                  }

                }
                else
                {
                    interfaceCell_[cellI]=false;
                }


            }
            else
            {

                normal_[cellI] = vector::zero;
                centre_[cellI] = vector::zero;
                interfaceCell_[cellI]=false;
            }
         }
         else
         {
            normal_[cellI] = vector::zero;
            centre_[cellI] = vector::zero;
            interfaceCell_[cellI]=false;
         }


    }



    if (writeVTK_ & mesh_.time().writeTime())
    {
        std::ostringstream os ;
        os << "isoFaces_" << int(mesh_.time().timeIndex());
        isoFacesToFile(facePts, os.str() , "isoFaces");
    }


}
