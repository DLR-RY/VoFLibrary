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

Application
    reconstructError

Description
    Reconstructs a VoF Field with a local centre and normal from a given
    Volume of fluid Field

Author:
    Henning Scheufler, DLR, all rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "OFstream.H"

#include "reconstructionSchemes.H"
#include "implicitFunctions.H"
#include "cutCellIso.H"
#include "reconstructionError.H"
#include "isoCutCell.H"




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


void setAlpha(const fvMesh& mesh, const dictionary& initAlphaFieldDict, volScalarField& alpha1)
{


    Foam::autoPtr<Foam::implicitFunctions> func =  implicitFunctions::New
    (
           initAlphaFieldDict.get<word>("function"),
           initAlphaFieldDict
    );

    scalarField f(mesh.nPoints(),0.0);

    forAll(f,pI)
    {
        f[pI] =  func->value(mesh.points()[pI]);
    }


    cutCellIso cutCell(mesh,f); // changes f
//     isoCutCell cutCell(mesh,f); // changes f


    forAll(alpha1,cellI)
    {
        label cellStatus = cutCell.calcSubCell(cellI,0.0);

        if(cellStatus == -1)
        {
            alpha1[cellI] = 1;
        }
        else if(cellStatus == 1)
        {
            alpha1[cellI] = 0;
        }
        else if(cellStatus == 0)
        {
            //alpha1[cellI]= max(min(cutCell.VolumeOfFluid(),1),0);
            alpha1[cellI]= max(min(cutCell.VolumeOfFluid(),1),0);
        }

    }

    alpha1.correctBoundaryConditions();

}

int main(int argc, char *argv[])
{
     #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
//    #include "createMesh.H"
     #include "createNamedMesh.H"


    IOdictionary initAlphaFieldDict
    (
        IOobject
        (
            "initAlphaFieldDict",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );



// init
    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    #include "createPhi.H"

    Info<< "Reading field alpha1\n" << endl;
    volScalarField alpha1
    (
        IOobject
        (
            "alpha1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    reconstructionError recErr( mesh,mesh,initAlphaFieldDict);


    #include "createTimeControls.H"
    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    //Setting velocity field and face fluxes for next time step
      Info << "Creating evaporationCurveModel\n" << endl;
      IOdictionary fvSolutionDict
      (
          IOobject
          (
              "fvSolution",
              alpha1.time().system(),
              alpha1.db(),
              IOobject::MUST_READ_IF_MODIFIED,
              IOobject::NO_WRITE
          )
      );

    runTime++;

    vector centre = initAlphaFieldDict.get<vector>("centre");
//    scalar radius = readScalar(initAlphaFieldDict.lookup("radius"));

    Random rndCentre(1234567);

    autoPtr<reconstructionSchemes> surf =  reconstructionSchemes::New(alpha1,phi,U,fvSolutionDict);

    label nIter = 100;

    word functionType = initAlphaFieldDict.get<word>("function");
    bool twoDim = (functionType == "disc");
    Info << "twoDim = " << twoDim << endl;
    Info << "functionType = " << functionType << endl;

    scalar recTime = 0;
    vector centreMin = centre - 0.1*centre;
    vector centreMax = centre + 0.1*centre;

    for(int iteration = 0;iteration < nIter;iteration++)
    {
        // 3D
        vector centrePos = rndCentre.globalPosition < vector > (centreMin,centreMax); //vector::zero;
//        if(twoDim)
//        {
//            centrePos = rndCentre.globalPosition < vector > (centreMin,centreMax);
//        }
//        else
//        {
//            centrePos = rndCentre.globalPosition < vector > (centreMin,centreMax);
//        }


        // 2D
//
        initAlphaFieldDict.set<vector>("centre",centrePos);

        centre = initAlphaFieldDict.lookup("centre");

//        Info << "centre " <<  centre << endl;


        setAlpha(mesh,initAlphaFieldDict,alpha1);
        mesh.time().cpuTimeIncrement();


        surf->reconstruct();

        recTime += mesh.time().cpuTimeIncrement();


        recErr.calcError(initAlphaFieldDict,surf->centre(),surf->normal(),false);




    }

    recErr.write(recTime/nIter,scalar(0));


    runTime.write();

	
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //