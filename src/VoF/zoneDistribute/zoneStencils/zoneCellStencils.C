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

#include "zoneCellStencils.H"
#include "syncTools.H"
#include "dummyTransform.H"
#include "emptyPolyPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(zoneCellStencils, 0);
}


Foam::autoPtr<Foam::indirectPrimitivePatch>
Foam::zoneCellStencils::nonEmptyFacesPatch() const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    label nNonEmpty = 0;

    for (const polyPatch& pp : patches)
    {
        if (!isA<emptyPolyPatch>(pp))
        {
            nNonEmpty += pp.size();
        }
    }
    labelList nonEmptyFaces(nNonEmpty);
    nNonEmpty = 0;

    for (const polyPatch& pp : patches)
    {
        if (!isA<emptyPolyPatch>(pp))
        {
            label facei = pp.start();

            forAll(pp, i)
            {
                nonEmptyFaces[nNonEmpty++] = facei++;
            }
        }
    }

    return autoPtr<indirectPrimitivePatch>::New
    (
        IndirectList<face>
        (
            mesh_.faces(),
            nonEmptyFaces
        ),
        mesh_.points()
    );
}

Foam::autoPtr<Foam::indirectPrimitivePatch>
Foam::zoneCellStencils::allCoupledFacesPatch() const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    label nCoupled = 0;

    for (const polyPatch& pp : patches)
    {
        if (pp.coupled())
        {
            nCoupled += pp.size();
        }
    }
    labelList coupledFaces(nCoupled);
    nCoupled = 0;

    for (const polyPatch& pp : patches)
    {
        if (pp.coupled())
        {
            label facei = pp.start();

            forAll(pp, i)
            {
                coupledFaces[nCoupled++] = facei++;
            }
        }
    }

    return autoPtr<indirectPrimitivePatch>::New
    (
        IndirectList<face>
        (
            mesh_.faces(),
            coupledFaces
        ),
        mesh_.points()
    );
}

void Foam::zoneCellStencils::validBoundaryFaces(boolList& isValidBFace) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    isValidBFace.setSize(mesh_.nFaces()-mesh_.nInternalFaces());

    isValidBFace = true;

    for (const polyPatch& pp : patches)
    {
        if (pp.coupled() || isA<emptyPolyPatch>(pp))
        {
            label bFacei = pp.start()-mesh_.nInternalFaces();
            forAll(pp, i)
            {
                isValidBFace[bFacei++] = false;
            }
        }
    }
}

void Foam::zoneCellStencils::unionEqOp::operator()
(
    labelList& x,
    const labelList& y
) const
{
    if (y.size())
    {
        if (x.empty())
        {
            x = y;
        }
        else
        {
            labelHashSet set(x);
            set.insert(y);
            x = set.toc();
        }
    }
}

void Foam::zoneCellStencils::merge
(
    const label globalI,
    const labelList& pGlobals,
    labelList& cCells
)
{
    labelHashSet set;
    for (const label celli : cCells)
    {
        if (celli != globalI)
        {
            set.insert(celli);
        }
    }

    for (const label celli : pGlobals)
    {
        if (celli != globalI)
        {
            set.insert(celli);
        }
    }

    cCells.setSize(set.size()+1);
    label n = 0;
    cCells[n++] = globalI;

    for (const label seti : set)
    {
        cCells[n++] = seti;
    }
}

void Foam::zoneCellStencils::insertFaceCells
(
    const label exclude0,
    const label exclude1,
    const boolList& isValidBFace,
    const labelList& faceLabels,
    labelHashSet& globals
) const
{
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    forAll(faceLabels, i)
    {
        label facei = faceLabels[i];

        label globalOwn = globalNumbering().toGlobal(own[facei]);
        if (globalOwn != exclude0 && globalOwn != exclude1)
        {
            globals.insert(globalOwn);
        }

        if (mesh_.isInternalFace(facei))
        {
            label globalNei = globalNumbering().toGlobal(nei[facei]);
            if (globalNei != exclude0 && globalNei != exclude1)
            {
                globals.insert(globalNei);
            }
        }
        else
        {
            label bFacei = facei-mesh_.nInternalFaces();

            if (isValidBFace[bFacei])
            {
                label globalI = globalNumbering().toGlobal
                (
                    mesh_.nCells()
                + bFacei
                );

                if (globalI != exclude0 && globalI != exclude1)
                {
                    globals.insert(globalI);
                }
            }
        }
    }
}


Foam::labelList Foam::zoneCellStencils::calcFaceCells
(
    const boolList& isValidBFace,
    const labelList& faceLabels,
    labelHashSet& globals
) const
{
    globals.clear();

    insertFaceCells
    (
        -1,
        -1,
        isValidBFace,
        faceLabels,
        globals
    );

    return globals.toc();
}




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneCellStencils::zoneCellStencils(const fvMesh& mesh)
:
    MeshObject<fvMesh, Foam::UpdateableMeshObject, zoneCellStencils>(mesh),
    labelListList(mesh.nCells()),
    mesh_(mesh),
    needComm_(),
    globalNumbering_(mesh_.nCells()+mesh_.nFaces()-mesh_.nInternalFaces())
{

}


void Foam::zoneCellStencils::updateMesh(const mapPolyMesh& mpm)
{
    if(mesh_.topoChanging())
    { 
        globalNumbering_ = globalIndex
        (
            mesh_.nCells()
            + mesh_.nFaces()
            - mesh_.nInternalFaces()
        );
        labelListList& stencil = *this;
        stencil = labelListList(mesh_.nCells());
        needComm_.clear();
    } 
}


// ************************************************************************* //
