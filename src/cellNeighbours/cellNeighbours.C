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

#include "cellNeighbours.H"
#include "dummyTransform.H"
#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"
#include "syncTools.H"
#include "wedgePolyPatch.H"

#include "globalPoints.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cellNeighbours::syncInterfacePoints(const boolList& nextToInterface)
{
    const labelListList& cPoints = mesh_.pointCells();
    const labelListList& pCells = mesh_.cellPoints();

    // includes non processor boundary faces
    syncInterfacePointMap_.clearStorage();

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    forAll(pbm, i) // mark on point patches
    {
        const polyPatch& pp = pbm[i];

        // exclude empty and wedgepatches
        if (!isA<emptyPolyPatch>(pp) && !isA<wedgePolyPatch>(pp))
        {
            const labelList& bPoints = pp.meshPoints();
            forAll(bPoints, j)
            {
                const label pI = bPoints[j];
                bool hasInterface = false;
                forAll(mesh_.pointCells()[pI], k)
                {
                    const label cellI = cPoints[pI][k];
                    if (nextToInterface[cellI])
                    {
                        hasInterface = true;
                        break;
                    }
                }
                if (hasInterface)
                {
                    syncInterfacePointMap_.insert(pI, hasInterface);
                }
            }
        }
    }

    // sync syncInterfacePointMap_ over processors
    syncTools::syncPointMap(mesh_, syncInterfacePointMap_, orEqOp<bool>());
}


void Foam::cellNeighbours::patchIFace
(
    const label faceI,
    label& patchNumber,
    label& patchFaceNumber
) const
{
    if (mesh_.isInternalFace(faceI))
    {
        patchNumber = -1;
    }
    else
    {
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

        // Boundary face. Find out which face of which patch
        const label patchI = pbm.patchID()[faceI - mesh_.nInternalFaces()];

        if (patchI < 0 || patchI >= pbm.size())
        {
            FatalErrorInFunction << "Cannot find patch for face " << faceI
                                 << abort(FatalError);
        }

        // Handle empty patches

        const polyPatch& pp = pbm[patchI];
        if (isA<emptyPolyPatch>(pp) || pp.empty() || isA<wedgePolyPatch>(pp) ||
            isA<processorPolyPatch>(pp))
        {
            patchNumber = -1;
        }
        else
        {
            patchNumber = patchI;
        }

        const label patchFaceI = pp.whichFace(faceI);
        patchFaceNumber = patchFaceI;
    }
}

void Foam::cellNeighbours::syncGlobalIdxMap(const boolList& nextToInterface)
{

    // sync syncInterfacePointMap_ on every processor
    // includes boundary faces
    syncInterfacePoints(nextToInterface);


    const labelListList& cPoints = mesh_.pointCells();
    const labelListList& fPoints = mesh_.pointFaces();

    // includes non processor boundary faces
    syncNeiCellMap_.clearStorage();

    DynamicList<label> IdxList(16); // 16 pointsNeibours

    // Insert all cells adjacent to the point in stencil
    forAllConstIters(syncInterfacePointMap_, iter)
    {
        IdxList.clear();
        const label pI = iter.key();
        forAll(cPoints[pI], j)
        {
            const label& pCellI = cPoints[pI][j];
            IdxList.append(globalNumbering_.toGlobal(pCellI));
        }
        // insert non wedge empty proc patches in stencil
        // if(!isInternalPoint(pI)) // check if it is internal point
        // {
        forAll(fPoints[pI], fi)
        {
            const label& faceI = fPoints[pI][fi];
            if (!mesh_.isInternalFace(faceI))
            {
                label bFacei = faceI - mesh_.nInternalFaces();
                label gIndex = globalNumbering_.toGlobal(mesh_.nCells() + bFacei);
                label patchI = -1;
                label i = -1;

                patchIFace(faceI, patchI, i);
                if (patchI != -1)
                {
                    IdxList.append(gIndex);//globalNumbering_.toGlobal(gIndex));
                }
            }
        }
        // }

        syncNeiCellMap_.insert(pI, IdxList);
    }

    // syncNeiCellMap processor values are appended to the list
    syncTools::syncPointMap(mesh_, syncNeiCellMap_, appendOp<label>());
}

void Foam::cellNeighbours::markPointPatches()
{
    onProcPatch_ = false;

    // set all points on proc patch to true
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    forAll(pbm, i)
    {
        const polyPatch& pp = pbm[i];

        if (isA<processorPolyPatch>(pp))
        {
            const labelList& bPoints = pp.meshPoints();
            forAll(bPoints, j)
            {
                const label pI = bPoints[j];
                onProcPatch_[pI] = true;
            }
        }
    }
}
void Foam::cellNeighbours::appendUniqueProcs
(
    const label& pI,
    DynamicList<label>& procs
)
{
    // pI the pointlabel has to be on the processors
    if (procPointsMap_.found(pI)) // check if on proc patch
    {
        const label procListI = procPointsMap_[pI];
        forAll(procPoints_[procListI], i)
        {
            const labelPair procPair = procPoints_[procListI][i];
            if (Pstream::myProcNo() != procPair.second())
            {
                // check if already present
                if (procs.size() != 0)
                {
                    bool exist = false;
                    forAll(procs, j)
                    {
                        if (procPair.second() == procs[j])
                        {
                            exist = true;
                        }
                    }

                    if (!exist)
                    {
                        procs.append(procPair.second());
                    }
                }
                else
                {
                    procs.append(procPair.second());
                }
            }
        }
    }
}

void Foam::cellNeighbours::calcBoundaryProcAddressing(const label& cellI)
{
    bool onProc = false;
    const labelListList& pCells = mesh_.cellPoints();
    const labelListList& cPoints = mesh_.pointCells();
    cellSet_.clearStorage();
    DynamicList<label> attachedProcessors(8);

    forAll(pCells[cellI], i)
    {
        const label& pI = pCells[cellI][i];
        if (onProcPatch_[pI]) // only processors
        {
            onProc = true;
            appendUniqueProcs(pI,attachedProcessors);
        }
    }

    if (onProc)
    {
        boundaryCellAddressing_.set
        (
            globalNumbering_.toGlobal(cellI),
            attachedProcessors
        ); // will  overwrite the entry
    }
}

void Foam::cellNeighbours::calcLocalCellStencil(const label& cellI)
{
    const labelListList& pCells = mesh_.cellPoints();
    const labelListList& cPoints = mesh_.pointCells();
    cellSet_.clearStorage();

    label globalCellI = globalNumbering_.toGlobal(cellI);
    cellSet_.insert(globalCellI);
    stencil_[cellI].append(globalCellI);

    // insert points map first
    bool onProc = false;

    forAll(pCells[cellI], i)
    {
        const label& pI = pCells[cellI][i];
        forAll(cPoints[pI], j)
        {

            if (syncNeiCellMap_.found(pI)) // only processors
            {
                if (onProcPatch_[pI]) // only processors
                {
                    onProc = true;
                }

                const labelList& index = syncNeiCellMap_[pI];
                if (index.size() != 0)
                {

                    forAll(index, i)
                    {
                        // is this slow cause
                        // every time we need to
                        // look
                        // may be it it faster to copy the labelList should
                        // only be 8 ish max
                        const label& gblIdx = index[i];
                        if (!cellSet_[gblIdx])
                        {
                            cellSet_.insert(gblIdx);
                            stencil_[cellI].append(gblIdx);

                        }
                    }
                }

            }
            else
            {
                // insert all adjacent cells
                const label& pCellI = cPoints[pI][j];
                label globalIdx = globalNumbering_.toGlobal(pCellI);
                if (!cellSet_[globalIdx])
                {
                    cellSet_.insert(globalIdx);
                    stencil_[cellI].append(globalIdx);
                }
            }
        }
    }

    if (onProc) // update boundaryCellAddressing_ if on processor
    {
        calcBoundaryProcAddressing(cellI);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellNeighbours::cellNeighbours(const fvMesh& mesh)
:
    mesh_(mesh),
    stencil_(mesh.nCells()),
    onProcPatch_(mesh.nPoints(), false),
    procPointsMap_(0),
    procPoints_(0),
    cellSet_(50),
    syncInterfacePointMap_(50),
    syncNeiCellMap_(mesh_.nFaces() - mesh_.nInternalFaces()), // is huge
    boundaryCellAddressing_(128), // will automatically update if too small
    globalNumbering_(mesh_.nCells() + mesh_.nFaces() - mesh_.nInternalFaces())
{

    globalPoints parallelPoints(mesh_, true, true);
    procPointsMap_ = parallelPoints.meshToProcPoint();
    procPoints_ = parallelPoints.procPoints();

    markPointPatches();

    forAll(stencil_, cellI)
    {
        stencil_[cellI].setCapacity(50);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * /5/

void Foam::cellNeighbours::updateStencil(const boolList& nextToInterface)
{
    const labelListList& cPoints = mesh_.pointCells();

    // sync BoundaryPoints
    // sets syncInterfacePointMap_
    // all on the boundary and nextToInterface are updated
    // syncInterfacePointMap_ in interface pointmap the globalindex are stored
    // including boundary faces without proc wedge empty patches
    syncGlobalIdxMap(nextToInterface);

    // the stencil is bigger than nextToInterface
    // loop over all marked boundary points
    // and check in which processor it is
    forAllConstIters(syncInterfacePointMap_, iter)
    {
        const label pI = iter.key();
        forAll(cPoints[pI], j)
        {
            const label& pCellI = cPoints[pI][j];
            if (!nextToInterface[pCellI])
            {
                calcBoundaryProcAddressing(pCellI);
            }
        }
    }

    // calcLocalCellStencil for every marked cell (nextToInterface)
    forAll(nextToInterface, cellI)
    {
        if (nextToInterface[cellI])
        {
            if (stencil_[cellI].size() == 0)
            {
                calcLocalCellStencil(cellI);
            }
        }
    }
}

// ************************************************************************* //
