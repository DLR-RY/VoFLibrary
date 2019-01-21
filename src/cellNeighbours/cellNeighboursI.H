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
#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"
#include "wedgePolyPatch.H"
#include "syncTools.H"

#include "Pair.H"
#include "globalPoints.H"

// ************************************************************************* //

template <typename Type>
void Foam::cellNeighbours::faceValue
(
    const GeometricField<Type, fvPatchField, volMesh>& f,
    const label facei,
    label& patchNumber,
    label& patchFaceNumber
) const
{
    if (mesh_.isInternalFace(facei))
    {
        patchNumber = -1;
    }
    else
    {
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

        // Boundary face. Find out which face of which patch
        const label patchi = pbm.patchID()[facei - mesh_.nInternalFaces()];

        if (patchi < 0 || patchi >= pbm.size())
        {
            FatalErrorInFunction << "Cannot find patch for face " << facei
                                 << abort(FatalError);
        }

        // Handle empty patches

        const polyPatch& pp = pbm[patchi];
        if (isA<emptyPolyPatch>(pp) || pp.empty() || isA<wedgePolyPatch>(pp) ||
            isA<processorPolyPatch>(pp))
        {
            patchNumber = -1;
        }
        else
        {
            patchNumber = patchi;
        }

        const label patchFacei = pp.whichFace(facei);
        patchFaceNumber = patchFacei;
        // return f.boundaryField()[patchi][patchFacei];
    }
}

template <typename Type>
const Foam::List<Foam::DynamicList<Foam::label> > &
Foam::cellNeighbours::getDataFromOtherProcs
(
    const GeometricField<Type, fvPatchField, volMesh>& f,
    const boolList& nextToInterface,
    Map<Type>& map
)
{

    const labelListList& cPoints = mesh_.pointCells();
    const labelListList& fPoints = mesh_.pointFaces();

    updateStencil(nextToInterface); // update the stencil

    List<Map<Type> > send(Pstream::nProcs());

    // cannot be on empty or
    // wedge since they are
    // ignored in
    // syncGlobalIdxMap
    forAllConstIters(syncInterfacePointMap_, iter)
    {
        const label pI = iter.key(); // fairly sure

        if (onProcPatch_[pI])
        {
            forAll(cPoints[pI], ci)
            {
                const label& cellI = cPoints[pI][ci]; //
                label gblIdx = globalNumbering_.toGlobal(cellI);
                {
                    // will crash if value not supplied
                    forAll(boundaryCellAddressing_[gblIdx], j)
                    {
                        const label procI = boundaryCellAddressing_[gblIdx][j];
                        send[procI].insert(gblIdx, f[cellI]);

                        forAll(fPoints[pI], fi)
                        {
                            // global face find
                            // boundary values and get
                            // value
                            const label& faceI = fPoints[pI][fi];
                            if (!mesh_.isInternalFace(faceI))
                            {
                                label bFacei = faceI - mesh_.nInternalFaces();
                                label gIndex = globalNumbering_.toGlobal(mesh_.nCells() + bFacei);
                                label patchi = -1;
                                label i = -1;

                                //faceValue<Type>(f, faceI, patchi, i);
                                patchIFace(faceI, patchi, i);
                                if (patchi != -1)
                                {
                                    Type val = f.boundaryField()[patchi][i];
                                    send[procI].insert(gIndex, val);
                                }
                            }
                        }
                    }
                }
            }


        }
        else // if on non processor BC append faces
        {
            forAll(fPoints[pI], fi)
            {
                // global face find
                // boundary values and get
                // value
                const label& faceI = fPoints[pI][fi];
                if (!mesh_.isInternalFace(faceI))
                {
                    label bFacei = faceI - mesh_.nInternalFaces();
                    label gIndex = globalNumbering_.toGlobal(mesh_.nCells() + bFacei);
                    label patchi = -1;
                    label i = -1;

                    patchIFace(faceI, patchi, i);
                    if (patchi != -1)
                    {
                        Type val = f.boundaryField()[patchi][i];
                        map.insert(gIndex, val);
                    }
                }
            }
        }
    }

    if (Pstream::parRun())
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Stream data into buffer
        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            if (domain != Pstream::myProcNo()) // && send[domain].size())
            {
                // Put data into send buffer
                UOPstream toDomain(domain, pBufs);

                toDomain << send[domain];
            }
        }

        // wait until everything is written.
        pBufs.finishedSends();

        Map<Type> neiMap;

        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            neiMap.clearStorage();

            if (domain != Pstream::myProcNo())
            {
                // get data from send buffer
                UIPstream fromDomain(domain, pBufs);

                fromDomain >> neiMap;
                // append neiMap to boundaryAndProcData
                forAllConstIters(neiMap, iter)
                {
                    map.insert(iter.key(), iter());
                }
            }
        }
    }

    return stencil_;
}

template <typename Type> void
Foam::cellNeighbours::getPointDataFromOtherProcs
(
    const GeometricField<Type, fvPatchField, volMesh>& f,
    const boolList& interfacePoints,
    Map<List<Type> >& map
)
{

    const labelListList& cPoints = mesh_.pointCells();
    const labelListList& fPoints = mesh_.pointFaces();

    // cannot be on empty or
    // wedge since they are
    // ignored in
    // syncGlobalIdxMap
    DynamicList<Type> valList(16); // 16 pointsNeibours

    // Insert all cells adjacent to the point in stencil
    forAllConstIters(syncInterfacePointMap_, iter)
    {
        valList.clear();
        const label pI = iter.key();
        forAll(cPoints[pI], j)
        {
            const label& pCellI = cPoints[pI][j];
            valList.append(f[pCellI]);
        }
        // insert non wedge empty proc patches in stencil
            forAll(fPoints[pI], fi)
            {
                const label& faceI = fPoints[pI][fi];
                if (!mesh_.isInternalFace(faceI))
                {
                    label bFacei = faceI - mesh_.nInternalFaces();
                    label patchI = -1;
                    label i = -1;

                    patchIFace(faceI, patchI, i);
                    if (patchI != -1)
                    {
                        Type val = f.boundaryField()[patchI][i];
                        valList.append(val);
                    }
                }
            }


        map.insert(pI, valList);
    }

    // syncNeiCellMap processor values are appended to the list
    syncTools::syncPointMap(mesh_, map, appendOp<Type>());

}

template <typename Type>
const Foam::List<Foam::DynamicList<Foam::label> >&
Foam::cellNeighbours::getDataFromOtherProcs
(
    const PtrList<GeometricField<Type, fvPatchField, volMesh> >& f,
    const boolList& nextToInterface,
    Map<List<Type> >& map
)
{
    // check if PtrList and Map list have the same size not clue how todo that
    //  If(f.size() != map)

    const labelListList& cPoints = mesh_.pointCells();
    const labelListList& fPoints = mesh_.pointFaces();

    updateStencil(nextToInterface); // update the stencil

    List<Map<List<Type > > > send(Pstream::nProcs());

    List<Type> appendVal(f.size());

    // cannot be on empty or
    // wedge since they are
    // ignored in
    // syncGlobalIdxMap
    forAllConstIters(syncInterfacePointMap_, iter)
    {
        const label pI = iter.key(); // fairly sure

        // if on proc append cells
        if (onProcPatch_[pI])
        {
            forAll(cPoints[pI], ci)
            {
                const label& cellI = cPoints[pI][ci]; //
                label gblIdx = globalNumbering_.toGlobal(cellI);
                {
                    // crashes if index not found
                    forAll(boundaryCellAddressing_[gblIdx], j)
                    {
                        const label procI = boundaryCellAddressing_[gblIdx][j];
                        forAll(f, i)
                        {
                            appendVal[i] = f[i][cellI];
                        }

                        send[procI].insert(gblIdx, appendVal);

                        forAll(fPoints[pI], fi)
                        {
                            // global face find
                            // boundary values and get
                            // value

                            const label& faceI = fPoints[pI][fi];
                            if (!mesh_.isInternalFace(faceI))
                            {
                                label bFacei = faceI - mesh_.nInternalFaces();
                                label gIndex = globalNumbering_.toGlobal(mesh_.nCells() + bFacei);
                                label patchi = -1;
                                label idx = -1;
                                forAll(f, i)
                                {
                                    //faceValue<Type>(f[i], faceI, patchi, idx);
                                    patchIFace(faceI, patchi, idx);
                                    if (patchi != -1)
                                    {

                                        Type val = f[i].boundaryField()[patchi][idx];
                                        appendVal[i] = val;
                                    }
                                }

                                if (patchi != -1) // should be the same as above
                                {
                                    // append
                                    send[procI].insert(gIndex, appendVal);
                                }
                            }
                        }
                    }
                }
            }


        }
        else // if on nonproc bc append faces
        {
            forAll(fPoints[pI], fi)
            {
                // global face find
                // boundary values and get
                // value

                const label& faceI = fPoints[pI][fi];
                if (!mesh_.isInternalFace(faceI))
                {
                    label bFacei = faceI - mesh_.nInternalFaces();
                    label gIndex = globalNumbering_.toGlobal(mesh_.nCells() + bFacei);
                    label patchi = -1;
                    label idx = -1;
                    forAll(f, i)
                    {
                        patchIFace(faceI, patchi, idx);
                        if (patchi != -1)
                        {
                            Type val = f[i].boundaryField()[patchi][idx];
                            appendVal[i] = val;
                        }
                    }
                    if (patchi != -1) // should be the same as above
                    {
                        map.insert(gIndex, appendVal);
                    }
                }
            }
        }
    }



    if (Pstream::parRun())
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Stream data into buffer
        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {

            if (domain != Pstream::myProcNo()) // && send[domain].size())
            {
                // Put data into send buffer
                UOPstream toDomain(domain, pBufs);

                toDomain << send[domain];
            }
        }

        // wait until everything is written.
        pBufs.finishedSends();

        Map<List<Type> > neiMap;

        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            neiMap.clearStorage();

            if (domain != Pstream::myProcNo())
            {
                // get data from send buffer
                UIPstream fromDomain(domain, pBufs);

                fromDomain >> neiMap;

                forAllConstIters(neiMap, iter)
                {
                    map.insert(iter.key(), iter());
                }
            }
        }
    }



    return stencil_;
}

// ************************************************************************* //
