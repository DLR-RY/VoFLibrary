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

#include "reconstructionError.H"
#include "addToRunTimeSelectionTable.H"

#include "reconstructionSchemes.H"
#include "implicitFunctions.H"



// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::reconstructionError::writeFileHeader(Ostream& os) const
{
    // Add headers to output data
    writeHeader(os, "reconstructionError ");
    writeCommented(os, "Time");
    writeTabbed(os, "reconstructionError");
    os  << endl;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reconstructionError::reconstructionError
(
    const fvMesh& mesh,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    writeFile(obr, "reconstructionError", "error", dict),
    mesh_(mesh),
    LNormalInf_(0),
    LNormal1_(0),
    LCentreInf_(0),
    LCentre1_(0),
    avgLNormal_(100),
    avgLCentre_(100)
{


    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reconstructionError::~reconstructionError()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::reconstructionError::calcError(const dictionary& dict,const volVectorField& centre, const volVectorField& normal, bool randomMode)
{

    if(randomMode)
    {
        avgLCentre_.clear();
        avgLNormal_.clear();
        LNormalInf_ = 0;
        LNormal1_ = 0;
    }

    const Foam::autoPtr<Foam::implicitFunctions>& func =  implicitFunctions::New
    (
           dict.lookup("function"),
           dict
    );

    DynamicList<scalar> LNormal(1000);
    DynamicList<scalar> LCentre(1000);

    forAll(normal,cellI)
    {
        if(mag(normal[cellI]) != 0)
        {
           vector centreCellI = centre[cellI];
           vector normalCellI = -normal[cellI];
           normalCellI /= mag(normalCellI);

           vector exactNormal = func->grad(centreCellI);
           exactNormal /= mag(exactNormal);

           LNormal.append(1-(normalCellI & exactNormal));
           LCentre.append(mag(func->distanceToSurfaces(centreCellI)));
        }
    }


    // normal errors
    scalar Lmax = max(LNormal);
    scalar L1 = average(LNormal);

    // centre errors
    scalar LCmax = max(LCentre);
    scalar LC1 = average(LCentre);

    avgLNormal_.append(L1);
    avgLCentre_.append(LC1);

    if(Lmax > LNormalInf_)
    {
        LNormalInf_ = Lmax;
    }

    if(LCmax > LCentreInf_)
    {
        LCentreInf_ = LCmax;
    }

    LNormal1_ = average(avgLNormal_);
    LCentre1_ = average(avgLCentre_);


}


void Foam::reconstructionError::write(const scalar k1, const scalar k2)
{
    scalar length = average(mag(mesh_.delta())()).value();
    scalar nCells = mesh_.nCells();
    scalar avgVol = average(mesh_.V()).value();


    if (Pstream::master())
    {
        writeTime(file());

        file()
            << token::TAB << LNormalInf_
            << token::TAB << LNormal1_
            << token::TAB << LCentreInf_
            << token::TAB << LCentre1_
            << token::TAB << k1
            << token::TAB << k2
            << token::TAB << length
            << token::TAB << nCells
            << token::TAB << avgVol
            << endl;
    }


}

void Foam::reconstructionError::write()
{
    scalar length = average(mag(mesh_.delta())()).value();
    scalar nCells = mesh_.nCells();
    scalar avgVol = average(mesh_.V()).value();


    if (Pstream::master())
    {
        writeTime(file());

        file()
            << token::TAB << LNormalInf_
            << token::TAB << LNormal1_
            << token::TAB << LCentreInf_
            << token::TAB << LCentre1_
            << token::TAB << length
            << token::TAB << nCells
            << token::TAB << avgVol
            << endl;
    }


}


// ************************************************************************* //
