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


#include "plicRDF.H"
#include "addToRunTimeSelectionTable.H"

#include "interpolationCellPoint.H"
#include "fvc.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reconstruction
{
    defineTypeNameAndDebug(plicRDF, 0);
    addToRunTimeSelectionTable(reconstructionSchemes,plicRDF, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reconstruction::plicRDF::plicRDF
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


    vof2IsoTol_(readScalar(modelDict().lookup("vof2IsoTol" ))),
    surfCellTol_(readScalar(modelDict().lookup("surfCellTol" ))),
    tol_(modelDict().lookupOrDefault("tol" ,1e-6)),
    relTol_(modelDict().lookupOrDefault("relTol" ,0.1)),
    iteration_(modelDict().lookupOrDefault("iterations" ,5)),
    interpolateNormal_(modelDict().lookupOrDefault("interpolateNormal" ,true)),
    heightFunc_(mesh_,  alpha1 ,surfCellTol_),
    sIterPLIC_(mesh_,surfCellTol_)

{
  writeVTK_ =  readBool( modelDict().lookup("writeVTK" ));
  heightFunc_.markCellsNearSurf();
  interfaceNormal_ = fvc::grad(alpha1_);

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
              interfaceCell_[cellI]=true;

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
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reconstruction::plicRDF::~plicRDF()
{}

// * * * * * * * * * * * * * * Protected Access Member Functions  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //
void Foam::reconstruction::plicRDF::reconstruct()
{

    if(interpolateNormal_)
    {
        Info << "interpolating normal" << endl;
        heightFunc_.markCellsNearSurf();

        heightFunc_.interpolateNormals(centre_,  normal_,U_,interfaceNormal_);
    }
    else
    {
        heightFunc_.markCellsNearSurf();
        heightFunc_.grad(alpha1_,interfaceNormal_);
    }

  
    DynamicList< List<point> > facePts;

    Map<bool> tooCoarse;

    for(int iter=0;iter<iteration_;iter++)
    {
        forAll(alpha1_,cellI)
        {
            if(sIterPLIC_.isASurfaceCell(alpha1_[cellI]))
            {
                if(!tooCoarse.found(cellI))
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
                          if (writeVTK_ & mesh_.time().writeTime() && iter==(iteration_-1))
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
        Map<scalar> residual;
        Map<scalar> avgAngle;

        //if(iter < (iteration_-1))
        {

            heightFunc_.calcHeightFunc(centre_,normal_,height_);
            heightFunc_.grad(height_,interfaceNormal_);
            heightFunc_.calcResidual(normal_,interfaceNormal_,residual,avgAngle);
        }


      label resCounter = 0;
      scalar avgRes = 0;
      scalar avgNormRes = 0;


      Map<scalar>::iterator resIter = residual.begin();
      Map<scalar>::iterator avgAngleIter = avgAngle.begin();

      while(resIter != residual.end())
      {


            if(avgAngleIter() > 0.52 && iter > 0) // 30 deg
            {
                tooCoarse.set(resIter.key(),true);
            }
            else
            {
                avgRes += resIter();
                scalar normRes = 0;
                scalar discreteError = 0.01*sqr(avgAngleIter());
                if(discreteError != 0)
                {
                    normRes= resIter()/max(discreteError,tol_);
                }
                else
                {
                    normRes= resIter()/tol_;
                }
                avgNormRes += normRes;
                resCounter++;

            }


            ++resIter;
            ++avgAngleIter;
      }

      reduce(avgRes,sumOp<scalar>());
      reduce(avgNormRes,sumOp<scalar>());
      reduce(resCounter,sumOp<label>());

      if(resCounter == 0) // avoid division  by zero and leave loop
      {
          resCounter = 1;
          avgRes = 0;
          avgNormRes = 0;
      }




      if(iter == 0)
      {
          Info << "intial residual absolute = " << avgRes/resCounter << endl;
          Info << "intial residual normalized = " << avgNormRes/resCounter << endl;
      }


      if(((avgNormRes/resCounter < relTol_ || avgRes/resCounter < tol_) && (iter >= 1 )) || iter + 1  == iteration_)
      {
          Info << "iterations = " << iter << endl;
          Info << "final residual absolute = " << avgRes/resCounter << endl;
          Info << "final residual normalized = " << avgNormRes/resCounter << endl;
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

