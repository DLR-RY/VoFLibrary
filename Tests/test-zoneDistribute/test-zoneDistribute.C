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
    test-zoneDistribute

Description
    test of zoneDistribute validated with mapDistribute
    
Author
    Henning Scheufler, DLR, all rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "centredCPCCellToCellStencilObject.H"
#include "zoneDistribute.H"

#include "SortableList.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const extendedCentredCellToCellStencil& stencil = centredCPCCellToCellStencilObject::New(mesh);

    const List<List<label> >& stencilAddr = stencil.stencil();

    List<List<vector> > stencilCentre(mesh.nCells());

    stencil.collectData
    (
       mesh.C(),
       stencilCentre
    );

    zoneDistribute exchangeFields(mesh);


    boolList interfaceCell(mesh.nCells(),true);
    exchangeFields.setUpCommforZone(interfaceCell);
    

    const labelListList& stencil_zoneDist = exchangeFields.getStencil();
    const globalIndex& gblNumbering = exchangeFields.globalNumbering();

    Map<vectorField> mapCC(exchangeFields.getFields(interfaceCell,mesh.C()));

    // compare stencils
    Pout << "size of the stencil match " << (stencilAddr.size() == stencil_zoneDist.size()) << endl;
    label stencilMatch = 0;

    forAll(stencilAddr,celli)
    {
        const vectorField& neiCC = mapCC[celli];
        bool foundAllLabel = true;
        if(neiCC.size() != stencilCentre[celli].size())
        {
            continue;
        }

        forAll(neiCC,i)
        {
            vector cc = neiCC[i];
            if(!stencilCentre[celli].found(cc))
            {
                foundAllLabel=false;
            }
        }
        if(foundAllLabel)
        {
            stencilMatch++;
        }

    }
    
    if(stencilMatch == mesh.nCells())
    {
        Pout << "all Values are identical "  << endl;
    }
    else
    {
        Pout << "values didnot match in : " << stencilMatch << " of "
             << mesh.nCells() << " cases" << endl;
    }


    return 0;
}


// ************************************************************************* //
