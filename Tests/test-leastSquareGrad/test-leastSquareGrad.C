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
    isoSurf

Description
    Uses isoCutter to create a volume fraction field from either a cylinder, 
    a sphere or a plane.

Author
    Henning Scheufler, DLR, all rights reserved.

\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "leastSquareGrad.H"

#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Vector<label> geomDim(1,1,-1);

    leastSquareGrad<scalar> lsGrad("polyDegree1",geomDim);

    List<vector> pos(6);
    List<scalar>  values(6);

    pos[0] = vector(0,0,0);
    pos[1] = vector(1,0,0);
    pos[2] = vector(-1,0,0);
    pos[3] = vector(0,1,0);
    pos[4] = vector(0,-1,0);
    //pos[5] = vector(1,1,0);

    values[0] = 0;
    values[1] = 1;
    values[2] = -1;
    values[3] = 1;
    values[4] = -1;
    //values[5] = 2;

    Map < List<vector>> posMap;
    Map < List<scalar>> valMap;

    posMap.insert(0,pos);
    posMap.insert(1,pos);
    posMap.insert(2,pos);

    valMap.insert(0,values);
    valMap.insert(1,values);
    valMap.insert(2,values);

    Info << lsGrad.grad(posMap,valMap) << endl;


    leastSquareGrad<vector> lsGradVec("polyDegree1",geomDim);

    List<vector>  valuesVec(6);

    valuesVec[0] = vector(0,0,0);
    valuesVec[1] = vector(1,0,0);
    valuesVec[2] = vector(-1,0,0);
    valuesVec[3] = vector(1,0,0);
    valuesVec[4] = vector(-1,0,0);

    Map < List<vector>> valMapVec;


    valMapVec.insert(0,valuesVec);
    valMapVec.insert(1,valuesVec);
    valMapVec.insert(2,valuesVec);

    Info << lsGradVec.grad(posMap,valMapVec) << endl;


    return 0;
}


// ************************************************************************* //
