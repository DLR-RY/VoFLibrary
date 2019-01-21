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


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reconstruction
{
    defineTypeNameAndDebug(gradAlpha, 0);
    addToRunTimeSelectionTable(reconstructionSchemes,gradAlpha, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reconstruction::gradAlpha::gradAlpha
(
        const volScalarField& alpha1,
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
    sIterPLIC_(mesh_,surfCellTol_),
    heightFunc_(mesh_,  alpha1 ,surfCellTol_)




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

    heightFunc_.markCellsNearSurf();
    heightFunc_.grad(alpha1_,interfaceNormal_);


    DynamicList< List<point> > facePts;

    forAll(alpha1_,cellI)
    {
        if(sIterPLIC_.isASurfaceCell(alpha1_[cellI]))
        {
               sIterPLIC_.vofCutCell
               (
                       cellI,
                       alpha1_[cellI],
                       vof2IsoTol_,
                       100,
                       interfaceNormal_[cellI]
               );


               if(sIterPLIC_.cellStatus() == 0)
               {

                   normal_[cellI] = sIterPLIC_.surfaceArea();
                   centre_[cellI] = sIterPLIC_.surfaceCentre();

                   if(mag(normal_[cellI]) != 0)
                   {
                        interfaceCell_[cellI]=true;
                        if (writeVTK_ & mesh_.time().writeTime())
                        {
                            facePts.append(sIterPLIC_.facePoints());
                        }
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

