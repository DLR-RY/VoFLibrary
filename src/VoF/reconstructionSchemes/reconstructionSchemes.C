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

#include "reconstructionSchemes.H"
#include "OFstream.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(reconstructionSchemes, 0);
    defineRunTimeSelectionTable(reconstructionSchemes, components);
}

void Foam::reconstructionSchemes::isoFacesToFile
(
    const DynamicList< List<point> >& faces,
    const word filNam,
    const word filDir
) const
{
    //Writing isofaces to vtk file for inspection in paraview

    mkDir(filDir);
    autoPtr<OFstream> vtkFilePtr;
    if (Pstream::parRun())
    {
        // Collect points from all the processors
        List<DynamicList<List<point>>> allProcFaces(Pstream::nProcs());
        allProcFaces[Pstream::myProcNo()] = faces;
        Pstream::gatherList(allProcFaces);

        if (Pstream::master())
        {
            Info << "Writing file: " << (filDir + "/" + filNam + ".vtk") << endl;
            vtkFilePtr.reset(new OFstream(filDir + "/" + filNam + ".vtk"));
            vtkFilePtr() << "# vtk DataFile Version 2.0" << endl;
            vtkFilePtr() << filNam << endl;
            vtkFilePtr() << "ASCII" << endl;
            vtkFilePtr() << "DATASET POLYDATA" << endl;
            //label nPoints(0);

            face f;
            label nPoints(0);
            label fSize = 0;
            forAll(allProcFaces, proci)
            {
                const DynamicList<List<point>>& procFaces = allProcFaces[proci];


                forAll(procFaces,fi)
                {
                    nPoints += procFaces[fi].size();
                    fSize++ ;
                }

            }

            vtkFilePtr() << "POINTS " << nPoints << " float" << endl;

            forAll(allProcFaces, proci)
            {
                const DynamicList<List<point>>& procFaces = allProcFaces[proci];

                forAll(procFaces,fi)
                {
                    List<point> pf = procFaces[fi];
                    forAll(pf,pi)
                    {
                        point p = pf[pi];
                        vtkFilePtr() << p[0] << " " << p[1] << " " << p[2] << endl;
                    }
                }

            }

            vtkFilePtr() << "POLYGONS " << fSize << " " << nPoints + fSize << endl;

            label np = 0;

            forAll(allProcFaces, proci)
            {
                const DynamicList<List<point>>& procFaces = allProcFaces[proci];


                forAll(procFaces,fi)
                {
                    nPoints = procFaces[fi].size();
                    vtkFilePtr() << nPoints;
                    for (label pi = np; pi < np + nPoints; pi++ )
                    {
                        vtkFilePtr() << " " << pi;
                    }
                    vtkFilePtr() << "" << endl;
                    np += nPoints;
                }

            }
        }
    }
    else
    {
        Info << "Writing file: " << (filDir + "/" + filNam + ".vtk") << endl;
        vtkFilePtr.reset(new OFstream(filDir + "/" + filNam + ".vtk"));
        vtkFilePtr() << "# vtk DataFile Version 2.0" << endl;
        vtkFilePtr() << filNam << endl;
        vtkFilePtr() << "ASCII" << endl;
        vtkFilePtr() << "DATASET POLYDATA" << endl;
        label nPoints(0);
        forAll(faces,fi)
        {
            nPoints += faces[fi].size();
        }

        vtkFilePtr() << "POINTS " << nPoints << " float" << endl;
        forAll(faces,fi)
        {
            List<point> pf = faces[fi];
            forAll(pf,pi)
            {
                point p = pf[pi];
                vtkFilePtr() << p[0] << " " << p[1] << " " << p[2] << endl;
            }
        }
        vtkFilePtr() << "POLYGONS " << faces.size() << " " << nPoints + faces.size() << endl;

        label np = 0;
        forAll(faces,fi)
        {
            nPoints = faces[fi].size();
            vtkFilePtr() << nPoints;
            for (label pi = np; pi < np + nPoints; pi++ )
            {
                vtkFilePtr() << " " << pi;
            }
            vtkFilePtr() << "" << endl;
            np += nPoints;
        }
    }
}

void Foam::reconstructionSchemes::isoFacesToFile
(
    const List<word> &fieldName,
    const List< DynamicList< double>>& values,
    const DynamicList< List<point> >& faces,
    const word filNam,
    const word filDir
) const
{
    //Writing isofaces to vtk file for inspection in paraview

    mkDir(filDir);
    autoPtr<OFstream> vtkFilePtr;
    Info << "Writing file: " << (filDir + "/" + filNam + ".vtk") << endl;
    vtkFilePtr.reset(new OFstream(filDir + "/" + filNam + ".vtk"));
    vtkFilePtr() << "# vtk DataFile Version 2.0" << endl;
    vtkFilePtr() << filNam << endl;
    vtkFilePtr() << "ASCII" << endl;
    vtkFilePtr() << "DATASET POLYDATA" << endl;

    label nPoints(0);
    forAll(faces,fi)
    {
        nPoints += faces[fi].size();
    }

    vtkFilePtr() << "POINTS " << nPoints << " float" << endl;
    forAll(faces,fi)
    {
        List<point> pf = faces[fi];
        forAll(pf,pi)
        {
            point p = pf[pi];
            vtkFilePtr() << p[0] << " " << p[1] << " " << p[2] << endl;
        }
    }
    vtkFilePtr() << "POLYGONS " << faces.size() << " " << nPoints + faces.size() << endl;

    label np = 0;
    forAll(faces,fi)
    {
        nPoints = faces[fi].size();
        vtkFilePtr() << nPoints;
        for (label pi = np; pi < np + nPoints; pi++ )
        {
            vtkFilePtr() << " " << pi;
        }
        vtkFilePtr() << "" << endl;
        np += nPoints;
    }


    vtkFilePtr() << "CELL_DATA " << faces.size()  << endl;
    vtkFilePtr() << "FIELD FieldData " << values.size()  << endl;
    forAll(values,i)
    {
        vtkFilePtr() << fieldName[i] << " 1 " <<  faces.size() << " double"  << endl;

        forAll(faces,fi)
        {
            vtkFilePtr() << values[i][fi] << endl;
        }
    }
}

bool Foam::reconstructionSchemes::alreadyReconstructed()
{
    const fvMesh& mesh = alpha1_.mesh();
    label& curTimeIndex = timeIndexAndIter_.first();
    label& curIter = timeIndexAndIter_.second();

    // rest timeIndex and curIter
    if(mesh.time().timeIndex() > curTimeIndex)
    {
        // maybe problematic with nOuterCorrectores
        // if(curIter >= 1)
        // {
        //     curIter = 0;
        //     return true;
        // }
        curTimeIndex = mesh.time().timeIndex();
        curIter = 0;
        return false;
    }
    
    // reconstruct always when subcycling
    if(mesh.time().subCycling() != 0)
    {
        return false;
    }
    
    curIter++;
    if(curIter > 1)
    {
        return true;
    }

    return false;

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reconstructionSchemes::reconstructionSchemes
(
        const word& type,
        volScalarField& alpha1,
        const surfaceScalarField& phi,
        const volVectorField& U,
        dictionary& dict
)
:
   IOdictionary
   (
        IOobject
        (
            "reconstructionScheme",
            alpha1.time().constant(),
            alpha1.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
  ),
  //reconstructionSchemesCoeffs_(subDict(type + "Coeffs")),
  reconstructionSchemesCoeffs_(dict),
  alpha1_(alpha1),
  phi_(phi),
  U_(U),
  normal_
  (
      IOobject
      (
          "recon::normal_",
          alpha1_.mesh().time().timeName(),
          alpha1_.mesh(),
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
      ),
      alpha1_.mesh(),
      dimensionedVector("0", dimArea, vector::zero),
      "calculated"
  ),
  centre_
  (
      IOobject
      (
          "recon::centre_",
          alpha1_.mesh().time().timeName(),
          alpha1_.mesh(),
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
      ),
      alpha1_.mesh(),
      dimensionedVector("0", dimLength, vector::zero),
      "calculated"
   ),
   interfaceCell_(alpha1_.mesh().nCells(),false),
   interfaceLabels_(0.2*alpha1_.mesh().nCells()),
   writeVTK_(false),
   timeIndexAndIter_(0,0)
{

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reconstructionSchemes::~reconstructionSchemes()
{}

// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * * * * * * * * //


const Foam::dictionary&
Foam::reconstructionSchemes::modelDict() const
{
    return reconstructionSchemesCoeffs_;
}

Foam::dictionary&
Foam::reconstructionSchemes::modelDict()
{
    return reconstructionSchemesCoeffs_;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //



// ************************************************************************* //
