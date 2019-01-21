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


#include "isoRDF.H"
#include "addToRunTimeSelectionTable.H"

#include "simpleMatrix.H"
#include "fvc.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reconstruction
{
    defineTypeNameAndDebug(isoRDF, 0);
    addToRunTimeSelectionTable(reconstructionSchemes,isoRDF, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reconstruction::isoRDF::isoRDF
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

    // Interpolation data

    vpi_(mesh_),
    ap_(mesh_.nPoints()),

    height_
    (
          IOobject
          (
              "height_",
              alpha1.mesh().time().timeName(),
              alpha1.mesh(),
              IOobject::NO_READ,
              IOobject::AUTO_WRITE
          ),
          sign(alpha1-0.5)*GREAT,
          "calculated"
    ),

    // Tolerances and solution controls
//    vof2IsoTol_(modelDict().lookupOrAddDefault<scalar>("vof2IsoTol", 1e-8,false,false)),
//    surfCellTol_(modelDict().lookupOrAddDefault<scalar>("surfCellTol", 1e-8,false,false))
    vof2IsoTol_(readScalar(modelDict().lookup("vof2IsoTol" ))),
    surfCellTol_(readScalar(modelDict().lookup("surfCellTol" ))),
    tol_(modelDict().lookupOrDefault("tol" ,1e-6)),
    iteration_(modelDict().lookupOrDefault("iterations" ,5)),
    heightFunc_(mesh_,  alpha1 ,surfCellTol_),
    sIterIso_(mesh_,ap_,surfCellTol_)




{
  writeVTK_ =  readBool( modelDict().lookup("writeVTK" ));
  Info << "iteration_ " << iteration_ << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reconstruction::isoRDF::~isoRDF()
{}

// * * * * * * * * * * * * * * Protected Access Member Functions  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //
void Foam::reconstruction::isoRDF::reconstruct()
{

    // Interpolating alpha1 cell centre values to mesh points (vertices)
    //ap_ = vpi_.interpolate(alpha1_);


    DynamicList< List<point> > facePts;

    heightFunc_.markCellsNearSurf();
    Map <vector> oldNormal;

    heightFunc_.interpolatePoints(alpha1_,ap_);

    for(int iter=0;iter<iteration_;iter++)
    {

    forAll(alpha1_,cellI)
    {
    //
        if(sIterIso_.isASurfaceCell(alpha1_[cellI]))
        {
            vector n(0,0,0);
            if(mag(normal_[cellI]) != 0)
            {
                n = normal_[cellI]/mag(normal_[cellI]);
            }


            oldNormal.set(cellI,n);



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
                  if (writeVTK_ & mesh_.time().writeTime()&& iter==(iteration_-1))
                  {
                          facePts.append(sIterIso_.facePoints());
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
         else
         {
            normal_[cellI] = vector::zero;
            centre_[cellI] = vector::zero;
            interfaceCell_[cellI]=false;
         }


    }

        normal_.correctBoundaryConditions();
        centre_.correctBoundaryConditions();


        heightFunc_.calcHeightFunc(centre_,normal_,height_);
        heightFunc_.interpolatePoints(height_,ap_);
        Map<scalar> residual;

        heightFunc_.calcResidual(normal_,oldNormal,residual);


        List<scalar> res(residual.size());
        label count = 0;
        forAllIter(Map<scalar>, residual, iter)
        {
            res[count++] = iter();
  //          Info << iter() << endl;
        }
        Info << "current residual absolute = " << gAverage(res) << endl;





        if(iter == 1)
        {
            Info << "intial residual absolute = " << gAverage(res) << endl;

        }

        if(( gAverage(res) < tol_ && (iter >= 1 )) || iter + 1  == iteration_)
        {
            Info << "iterations = " << iter << endl;
            Info << "final residual absolute = " << gAverage(res) << endl;
            break;
        }

    }


    if (writeVTK_ & mesh_.time().writeTime())
    {
        std::ostringstream os ;
        os << "isoFaces_" << int(mesh_.time().timeIndex());
        isoFacesToFile(facePts, os.str() , "isoFaces");
    }

}

