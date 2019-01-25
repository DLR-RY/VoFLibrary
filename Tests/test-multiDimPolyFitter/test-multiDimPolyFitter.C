/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of the IsoAdvector source code library, which is an 
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
    isoSurf

Description
    Uses isoCutter to create a volume fraction field from either a cylinder, 
    a sphere or a plane.

Author
    Henning Scheufler, DLR, all rights reserved.

\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "multiDimPolyFitter.H"

#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Vector<label> test(1,-1,-1);
    // test[0] = -1;
    // test[1] = 1;
    // test[2] = 1;
    autoPtr<multiDimPolyFunctions> polyFunc = multiDimPolyFunctions::New("polyDegree1",test);
    vector vec(1,2,3);

    Info << "polyFunc nTerms " << polyFunc->nTerms()  << endl;
    Info << "polyFunc coeffs " << polyFunc->coeffs()  << endl;
    Info << "polyFunc termValues " << polyFunc->termValues(vec)  << endl;

    multiDimPolyFitter < scalar > polyFitter("polyDegree1",test);

    List<vector> pos(3);
    List<scalar>  values(3);

    pos[0] = vector(0,0,0);
    pos[1] = vector(1,0,0);
    pos[2] = vector(-1,0,0);

    values[0] = 1;
    values[1] = 2;
    values[2] = 0;

    Info << "pos " << pos  << endl;
    Info << "values " << values  << endl;

    scalarField fitData = polyFitter.fitData
    (
        pos,
        values
    );

    Info << "fitData " << fitData  << endl;

    autoPtr<multiDimPolyFunctions > polyFuncDeg2 = multiDimPolyFunctions::New("polyDegree2",Vector<label>(1,1,-1));
    //vector vec(1,2,3);

    Info << "polyFunc nTerms " << polyFuncDeg2->nTerms()  << endl;
    // Info << "polyFunc coeffs " << polyFunc->coeffs()  << endl;
    Info << "polyFuncDeg2 termValues " << polyFuncDeg2->termValues(vector(1,1,0))  << endl;

    List<vector> pos2(6);
    List<scalar>  values2(6);

    pos2[0] = vector(0,0,0);
    pos2[1] = vector(1,0,0);
    pos2[2] = vector(-1,0,0);
    pos2[3] = vector(0,1,0);
    pos2[4] = vector(0,-1,0);
    pos2[5] = vector(1,1,0);

    values2[0] = 0;
    values2[1] = 1;
    values2[2] = 1;
    values2[3] = 1;
    values2[4] = 1;
    values2[5] = 2;

    multiDimPolyFitter < scalar > polyFitter2("polyDegree2",Vector<label>(1,1,-1));

    scalarField fitData2 = polyFitter2.fitData
    (
        pos2,
        values2
    );

    Info << "fitData " << fitData2  << endl;


    return 0;
}


// ************************************************************************* //
