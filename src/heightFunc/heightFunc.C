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

//#include "dummyTransform.H"
#include "emptyPolyPatch.H"
#include "heightFunc.H"
#include "processorPolyPatch.H"
#include "syncTools.H"
#include "wedgePolyPatch.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::heightFunc::loopStencil
(
    const DynamicList < label >& neiCells,
    const Map < List < vector> >& map,
    const globalIndex& globalIdx,
    const volVectorField& centre,
    const volVectorField& normal,
    const point& p,
    scalar& averageDist,
    scalar& avgWeight
)
{

    forAll(neiCells, i)
    {
        const label gIdx = neiCells[i];

        if (globalIdx.isLocal(gIdx))
        {
            const label localCI = globalIdx.toLocal(gIdx);
            if (localCI < mesh_.nCells()) // boundary values are not searched
            {
                vector n = -normal[localCI];

                if (mag(n) != 0)
                {
                    n /= mag(n);
                    vector distanceToIntSeg = (p - centre[localCI]);
                    scalar distToSurf = distanceToIntSeg & (n);
                    scalar weight = 0;

                    if (mag(distanceToIntSeg) != 0)
                    {
                        distanceToIntSeg /= mag(distanceToIntSeg);
                        weight = sqr(mag(distanceToIntSeg & n));
                    }
                    else // exactly on the center
                    {
                        weight = 1;
                    }
                    averageDist += distToSurf * weight;
                    avgWeight += weight;
                }
            }
        }
        else
        {

            const List<vector>& interface = map[gIdx]; // on procssor

            vector n = -interface[1];
            vector c = interface[0];
            if (mag(n) != 0)
            {
                n /= mag(n);
                vector distanceToIntSeg = (p - c);
                scalar distToSurf = distanceToIntSeg & (n);
                scalar weight = 0;

                if (mag(distanceToIntSeg) != 0)
                {
                    distanceToIntSeg /= mag(distanceToIntSeg);
                    weight = sqr(mag(distanceToIntSeg & n));
                }
                else // exactly on the center weight can be have a maximum of 1
                {
                    weight = 1;
                }
                averageDist += distToSurf * weight;
                avgWeight += weight;
            }
        }
    }
}

void Foam::heightFunc::loopStencilFindNearNeighbour
(
    const DynamicList<label>& neiCells,
    const Map<List<vector> >& map,
    const globalIndex& globalIdx,
    const volVectorField& centre,
    const volVectorField& normal,
    const point& p,
    vector& normalNearNei,
    DynamicList<vector>& interfaceNormals
)
{
    scalar smallDist = GREAT;

    forAll(neiCells, i)
    {

        const label gIdx = neiCells[i];

        if (globalIdx.isLocal(gIdx))
        {
            const label localCI = globalIdx.toLocal(gIdx);
            if (localCI < mesh_.nCells()) // boundary values are not weighted
            {
                vector n = -normal[localCI];

                if (mag(n) != 0)
                {
                    n /= mag(n);
                    vector distanceToIntSeg = (p - centre[localCI]);
                    if (mag(distanceToIntSeg) < smallDist)
                    {
                        smallDist = mag(distanceToIntSeg);
                        normalNearNei = n;

                    }
                    interfaceNormals.append(n);
                }
            }
        }
        else
        {
            // is on other proc
            // will crash if not found
            const List<vector>& interface = map[gIdx];

            vector n = -interface[1];
            vector c = interface[0];
            if (mag(n) != 0)
            {
                n /= mag(n);
                vector distanceToIntSeg = (p - c);

                if (mag(distanceToIntSeg) < smallDist)
                {
                    smallDist = mag(distanceToIntSeg);
                    normalNearNei = n;

                }
                interfaceNormals.append(n);
            }
        }
    }

}

void Foam::heightFunc::loopStencilInterpolateNormal
(
    const DynamicList<label>& neiCells,
    const Map<List<vector> >& map,
    const globalIndex& globalIdx,
    const volVectorField& centre,
    const volVectorField& normal,
    const point& p,
    vector& normalNearNei,
    DynamicList<vector>& interfaceNormals
)
{
    scalar smallDist = GREAT;
    scalar weight = 0;

    forAll(neiCells, i)
    {

        const label gIdx = neiCells[i];

        if (globalIdx.isLocal(gIdx))
        {
            const label localCI = globalIdx.toLocal(gIdx);
            if (localCI < mesh_.nCells()) // boundary values are not weighted
            {
                vector n = -normal[localCI];

                if (mag(n) != 0)
                {
                    n /= mag(n);
                    vector distanceToIntSeg = (tensor::I- n*n) & (p - centre[localCI]);// project vector in plane
                    normalNearNei += n /max(mag(distanceToIntSeg),SMALL);
                    weight += 1/max(mag(distanceToIntSeg),SMALL);

                    interfaceNormals.append(n);
                }
            }
        }
        else
        {
            // is on other proc
            // will crash if not found
            const List<vector>& interface = map[gIdx];

            vector n = -interface[1];
            vector c = interface[0];
            if (mag(n) != 0)
            {
                n /= mag(n);
                vector distanceToIntSeg = (tensor::I- n*n) & (p - c);// project vector in plane
                normalNearNei += n / max(mag(distanceToIntSeg),SMALL);
                weight += 1/max(mag(distanceToIntSeg),SMALL);

                interfaceNormals.append(n);
            }
        }
    }
    if(weight != 0)
    {
        normalNearNei /= weight;
    }

}

void Foam::heightFunc::loopStencilEstimatedNormalAndCentre
(
    const DynamicList<label>& neiCells,
    const Map<List<vector> >& map,
    const globalIndex& globalIdx,
    const volVectorField& centre,
    const volVectorField& normal,
    const point& p,
    vector& avgNormal,
    vector& interpolCentre
)
{
    scalar smallDist = GREAT;
    scalar avgWeight = 0;
    forAll(neiCells, i)
    {

        const label gIdx = neiCells[i];

        if (globalIdx.isLocal(gIdx))
        {
            const label localCI = globalIdx.toLocal(gIdx);
            if (localCI < mesh_.nCells()) // boundary values are not weighted
            {
                vector n = -normal[localCI];

                if (mag(n) != 0)
                {
                    n /= mag(n);
                    vector distanceToIntSeg = (p - centre[localCI]);
                    if(mag(distanceToIntSeg) < smallDist)
                    {
                        smallDist = mag(distanceToIntSeg);
                        interpolCentre = centre[localCI];
//                        avgNormal = n;

                    }
                    scalar weight = 0;

                    if (mag(distanceToIntSeg) != 0)
                    {
                        distanceToIntSeg /= mag(distanceToIntSeg);
                        weight = sqr(mag(distanceToIntSeg & n)); //sqr(mag(distanceToIntSeg & n));
                    }
                    else // exactly on the center
                    {
                        weight = 1;
                    }
                    avgNormal += n * weight;
                    avgWeight += weight;
                }
            }
        }
        else
        {
            // is on other proc
            // will crash if not found
            const List<vector>& interface = map[gIdx];

            vector n = -interface[1];
            vector c = interface[0];

            if (mag(n) != 0)
            {
                n /= mag(n);
                vector distanceToIntSeg = (p - c);
                if(mag(distanceToIntSeg) < smallDist)
                {
                    smallDist = mag(distanceToIntSeg);
                    interpolCentre = c;
//                    avgNormal = n;

                }
                scalar weight = 0;

                if (mag(distanceToIntSeg) != 0)
                {
                    distanceToIntSeg /= mag(distanceToIntSeg);
                    distanceToIntSeg /= mag(distanceToIntSeg);
                    weight = sqr(mag(distanceToIntSeg & n)); //sqr(mag(distanceToIntSeg & n));
                }
                else // exactly on the center
                {
                    weight = 1;
                }
                avgNormal += n * weight;
                avgWeight += weight;
            }
        }
    }

    if(avgWeight != 0)
    {
        avgNormal /= avgWeight;
    }
}

void Foam::heightFunc::loopStencilAdvectHeightAndEstimateNormal
(
    const DynamicList<label>& neiCells,
    const Map<List<vector> >& map,
    const globalIndex& globalIdx,
    const volVectorField& centre,
    const volVectorField& normal,
    const volVectorField& U,
    const point& p,
    scalar& height,
    vector& avgNormal,
    vector& interpolCentre
)
{
    scalar dt = mesh_.time().deltaTValue();
    scalar smallDist = GREAT;
    scalar avgWeight = 0;
    scalar avgAdvect = 0;
    forAll(neiCells, i)
    {

        const label gIdx = neiCells[i];

        if (globalIdx.isLocal(gIdx))
        {
            const label localCI = globalIdx.toLocal(gIdx);
            if (localCI < mesh_.nCells()) // boundary values are not weighted
            {
                vector n = -normal[localCI];

                if (mag(n) != 0)
                {
                    n /= mag(n);
                    vector distanceToIntSeg = (p - centre[localCI]);
//                    if(mag(distanceToIntSeg) < smallDist)
//                    {
//                        smallDist = mag(distanceToIntSeg);
//                        interpolCentre = centre[localCI];
//                    }

                    scalar weight = 0;

                    if (mag(distanceToIntSeg) != 0)
                    {
                        distanceToIntSeg /= mag(distanceToIntSeg);
                        weight = sqr(mag(distanceToIntSeg & n)); //sqr(mag(distanceToIntSeg & n));
                    }
                    else // exactly on the center
                    {
                        weight = 1;
                    }
                    avgNormal += n * weight;
                    avgWeight += weight;
                    avgAdvect += dt*(U[localCI] & n);
                }
            }
        }
        else
        {
            // is on other proc
            // will crash if not found
            const List<vector>& interface = map[gIdx];

            vector n = -interface[1];
            vector c = interface[0];

            if (mag(n) != 0)
            {
                n /= mag(n);
                vector distanceToIntSeg = (p - c);
//                if(mag(distanceToIntSeg) < smallDist)
//                {
//                    smallDist = mag(distanceToIntSeg);
//                    interpolCentre = c;
////                    avgNormal = n;

//                }
                scalar weight = 0;

                if (mag(distanceToIntSeg) != 0)
                {
                    distanceToIntSeg /= mag(distanceToIntSeg);
                    distanceToIntSeg /= mag(distanceToIntSeg);
                    weight = sqr(mag(distanceToIntSeg & n));
                }
                else // exactly on the center
                {
                    weight = 1;
                }
                avgNormal += n * weight;
                avgWeight += weight;
                //avgAdvect += dt*(U[localCI] & n);
            }
        }
    }

    if(avgWeight != 0)
    {
        avgNormal /= avgWeight;
        avgAdvect /= avgWeight;
    }

    height += avgAdvect;

    vector n = avgNormal/mag(avgNormal);

    interpolCentre = p - (n*height);

}


void Foam::heightFunc::markCellsNearSurf()
{

    const labelListList& pCells = mesh_.cellPoints();
    const labelListList& cPoints = mesh_.pointCells();
    interface_ = false;
    nextToInterface_ = false;
    interfacePoint_ = false;

    // do coupled face first
    Map<bool> syncMap;

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    forAll(pbm, i) // mark on point patches
    {
        const polyPatch& pp = pbm[i];

        if (!isA<emptyPolyPatch>(pp) && !isA<wedgePolyPatch>(pp))
        {
            const labelList& bPoints = pp.meshPoints();
            forAll(bPoints, j)
            {
                const label pI = bPoints[j];
                forAll(mesh_.pointCells()[pI], k)
                {
                    const label cellI = cPoints[pI][k];
                    if (isASurfaceCell(alpha1_[cellI]))
                    {
                        interfacePoint_[pI] = true;
                        break;
                    }
                }
                if (interfacePoint_[pI])
                {
                    syncMap.insert(pI, interfacePoint_[pI]);
                }
            }
        }
    }

    syncTools::syncPointMap(mesh_, syncMap, orEqOp<bool>());

    // loop Map and update nextToInterface_ and nextToInterface_
    forAllIter(Map<bool>, syncMap, iter)
    {
        interfacePoint_[iter.key()] = iter();
        forAll(cPoints[iter.key()],
               j) // loop over all cells attached to the point
        {
            const label& pCellI = cPoints[iter.key()][j];
            nextToInterface_[pCellI] = true;
        }
    }

    // update nextToInterface_ and interface_
    forAll(alpha1_, cellI)
    {
        if (isASurfaceCell(alpha1_[cellI]))
        {
            interface_[cellI] = true;
            forAll(pCells[cellI], i)
            {
                const label& pI = pCells[cellI][i];
                if (interfacePoint_[pI])
                {
                    // cell already marked do nothign
                }
                else
                {
                    interfacePoint_[pI] = true;
                    forAll(cPoints[pI],
                           j) // loop over all cells attached to the point
                    {
                        const label& pCellI = cPoints[pI][j];
                        nextToInterface_[pCellI] = true;
                    }
                }
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heightFunc::heightFunc
(
    const fvMesh& mesh,
    const volScalarField& alpha1,
    const scalar& tol
)
:
    mesh_(mesh),
    alpha1_(alpha1),
    surfCellTol_(tol),
    interface_(mesh.nCells(), false),
    nextToInterface_(mesh.nCells(), false),
    interfacePoint_(mesh.nPoints(), false),
    nei_(mesh)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::heightFunc::grad
(
    const volScalarField& height,
    volVectorField& surfnormals
)
{

    leastSquareCell lsGrad(mesh_, nei_.globalNumbering());

    // create map with height and cell centre information
    Map<scalar> mapHeight;
    Map<vector> mapCC;

    // fill maps and get stenicl
    const List<DynamicList<label> >& stencil =
        nei_.getDataFromOtherProcs < scalar > ( height, interface_, mapHeight);
    nei_.getDataFromOtherProcs < vector > ( mesh_.C(), interface_, mapCC);

    // loop cells and calculated the gradient in that cell
    forAll(surfnormals, cellI)
    {
        if (interface_[cellI])
        {
            vector normal = lsGrad.calcGrad(
                cellI, mesh_.C()[cellI], height, stencil[cellI], mapHeight, mapCC);

            if (mag(normal) != 0)
            {
                surfnormals[cellI] = normal;
            }
        }
    }
}

void Foam::heightFunc::grad
(
    const volScalarField& height,
    const volVectorField& centre,
    volVectorField& surfnormals
)
{

    leastSquareCell lsGrad(mesh_, nei_.globalNumbering());

    // create map with height and cell centre information
    Map<scalar> mapHeight;
    Map<vector> mapCC;

    // fill maps and get stenicl
    const List<DynamicList<label> >& stencil =
        nei_.getDataFromOtherProcs < scalar > ( height, interface_, mapHeight);
    nei_.getDataFromOtherProcs < vector > ( mesh_.C(), interface_, mapCC);

    // loop cells and calculated the gradient in that cell
    forAll(surfnormals, cellI)
    {
        if (interface_[cellI])
        {
            vector normal = lsGrad.calcGrad(
                cellI, centre[cellI], height, stencil[cellI], mapHeight, mapCC);

            if (mag(normal) != 0)
            {
                surfnormals[cellI] = normal;
            }
        }
    }
}

void Foam::heightFunc::calcResidual
(
        const volVectorField& normal, // old
        const volVectorField& interfaceNormal, // new
        Map<scalar>& normalResidual,
        Map<scalar>& avgAngle
)

{

    // create map with height and cell centre information
    Map<vector> mapNormal;
    const globalIndex& globalIdx = nei_.globalNumbering();


    // fill maps and get stenicl
    const List<DynamicList<label> >& stencil =
    nei_.getDataFromOtherProcs < vector > ( normal, interface_, mapNormal);

    // loop cells and calculated the gradient in that cell
    normalResidual.clear();

    label counter = 0;
    forAll(normal, cellI)
    {
        if (mag(normal[cellI]) != 0) //
        {
            scalar avgDiffNormal = 0;
            scalar maxDiffNormal = GREAT;
            scalar NofInterfaces= 0;
            const vector cellNormal = normal[cellI]/mag(normal[cellI]); // cache it
            forAll(stencil[cellI], i)
            {
                const label gIdx = stencil[cellI][i];
                if (globalIdx.isLocal(gIdx))
                {
                    // also includes boundary faces
                    const label localCellI = globalIdx.toLocal(gIdx);
                    if (localCellI <  mesh_.nCells())
                    {
                        if(mag(normal[localCellI]) != 0 && localCellI != cellI)
                        {
                            vector n = normal[localCellI]/mag(normal[localCellI]);
                            scalar cosAngle = max(min((cellNormal & n),1),-1);
//                            Info << "cosAngle " << cosAngle << endl;
                            avgDiffNormal += acos(cosAngle) * mag(normal[localCellI]);
                            NofInterfaces += mag(normal[localCellI]);
                            if(cosAngle < maxDiffNormal)
                            {
                                maxDiffNormal = cosAngle;
                            }
                        }
                    }
                }
                else
                {
                    vector n = mapNormal[gIdx];
                    if(mag(n) != 0)
                    {
                        n /= mag(n);
                    }
                    scalar cosAngle = max(min((cellNormal & n),1),-1);
                    if(cosAngle < maxDiffNormal)
                    {
                        maxDiffNormal = cosAngle;
                    }
                    //avgDiffNormal += (cellNormal & n);
                    //NofInterfaces++;
                    avgDiffNormal += acos(cosAngle) * mag(mapNormal[gIdx]);
                    NofInterfaces += mag(mapNormal[gIdx]) ;
                }
            }
            if(NofInterfaces != 0)
            {
                avgDiffNormal /= NofInterfaces;
            }
            else
            {
                avgDiffNormal = 0;
            }
            vector newCellNormal = -interfaceNormal[cellI]/mag(interfaceNormal[cellI]);
            scalar normalRes = (1 - (cellNormal & newCellNormal)); // max((1 - avgDiffNormal),1e-8);
            avgAngle.insert(cellI,avgDiffNormal);
            normalResidual.insert(cellI,normalRes);
        }
    }
}

void Foam::heightFunc::calcResidual
(
        const volVectorField& normal, // new
        Map<vector>& oldNormal, // old
        Map<scalar>& normalResidual
)

{

    // loop cells and calculated the gradient in that cell
    normalResidual.clear();

    label counter = 0;
    forAll(normal, cellI)
    {
        if (mag(normal[cellI]) != 0) //
        {
            const vector cellNormal = normal[cellI]/mag(normal[cellI]); // cache it
            vector oldCellNormal = oldNormal[cellI];
            if(oldCellNormal == vector::zero)
            {
                normalResidual.insert(cellI,1);
            }
            else
            {
                oldCellNormal /= mag(oldCellNormal);
                // normalResidual
                scalar normalRes = (1 - (cellNormal & oldCellNormal)); // max((1 - avgDiffNormal),1e-8);

                normalResidual.insert(cellI,normalRes);
            }
        }
    }
}

void Foam::heightFunc::interpolatePoints
(
    const volScalarField& height,
    scalarField& ap
)
{
    Info << "interpolate points" << endl;
    leastSquareCell lsInterpol(mesh_, nei_.globalNumbering());

    // create map with height and cell centre information
    Map <List <scalar> > mapHeight;
    Map <List <vector> > mapCC;


    nei_.getPointDataFromOtherProcs < scalar > ( height, interfacePoint_, mapHeight);
    nei_.getPointDataFromOtherProcs < vector > ( mesh_.C(), interfacePoint_, mapCC);

    // fill maps and get stenicl

    forAll(interfacePoint_,pI)
    {
        if(interfacePoint_[pI])
        {
                ap[pI] = lsInterpol.interpolate(pI, height,mapHeight,mapCC);
        }
    }

}


void Foam::heightFunc::calcHeightFunc
(
    const volVectorField& centre,
    const volVectorField& normal,
    volScalarField& height
)
{

    // create PtrList
    PtrList<volVectorField> interfaceData(2);
    // fill ptrlist
    interfaceData.set(0, centre);
    interfaceData.set(1, normal);

    // get filled by nei
    Map<List<vector> > map;

    const List<DynamicList<label> >& stencil =
        nei_.getDataFromOtherProcs<vector>(
            interfaceData, nextToInterface_, map); // interpolate a point value

    scalar smallDist = GREAT;
    scalar averageDist = 0;
    scalar avgWeight = 0;

    // get gobal addressing
    const globalIndex& globalIdx = nei_.globalNumbering();

    forAll(nextToInterface_, cellI)
    {
        if (nextToInterface_[cellI])
        {

            height[cellI] = sign(alpha1_[cellI] - 0.5) * GREAT;
            smallDist = GREAT;
            averageDist = 0;
            avgWeight = 0;

            loopStencil(stencil[cellI],
                        map,
                        globalIdx,
                        centre,
                        normal,
                        mesh_.C()[cellI],
                        averageDist,
                        avgWeight);

            if (avgWeight != 0)
            {
                height[cellI] = averageDist / avgWeight;
            }
        }
        else
        {
            height[cellI] = sign(alpha1_[cellI] - 0.5) * GREAT;
        }
    }

    forAll(height.boundaryField(), patchI)
    {
        if (height.boundaryField().types()[patchI] == "calculated")
        {
            const polyPatch pp = mesh_.boundaryMesh()[patchI];
            fvPatchScalarField& pHeight = height.boundaryFieldRef()[patchI];
            forAll(pHeight, i)
            {
                pHeight[i] =
                    sign(alpha1_.boundaryField()[patchI][i] - 0.5) * GREAT;
                const label& pCellI = pp.faceCells()[i];
                smallDist = GREAT;
                averageDist = 0;
                avgWeight = 0;

                if (nextToInterface_[pCellI])
                {

                    loopStencil
                    (
                        stencil[pCellI],
                        map,
                        globalIdx,
                        centre,
                        normal,
                        mesh_.C().boundaryField()[patchI][i],
                        averageDist,
                        avgWeight
                    );

                    if (avgWeight != 0)
                    {
                        pHeight[i] = averageDist / avgWeight;
                    }
                }
            }
        }
    }

    height.correctBoundaryConditions();
}

void Foam::heightFunc::interpolateNormals
(
    const volVectorField& centre,
    const volVectorField& normal,
    volVectorField& surfnormals
)
{

    leastSquareCell lsGrad(mesh_, nei_.globalNumbering());

    // create PtrList
    PtrList<volVectorField> interfaceData(2);
    // fill ptrlist
    interfaceData.set(0, centre);
    interfaceData.set(1, normal);

    // create maps
    Map<List<vector> > map;
    Map<scalar> mapAlpha;
    Map<vector> mapCC;

    // update maps and get stencil
    const List<DynamicList<label> >& stencil =
        nei_.getDataFromOtherProcs<scalar>(
            alpha1_, interface_, mapAlpha);
    // update maps
    nei_.getDataFromOtherProcs<vector>(mesh_.C(), interface_, mapCC);
    nei_.getDataFromOtherProcs<vector>(interfaceData, interface_, map);

    vector closestNormal = vector::zero;
    DynamicList<vector> foundNormals(10);

    const globalIndex& globalIdx = nei_.globalNumbering();

    forAll(alpha1_, cellI)
    {

        if (isASurfaceCell(alpha1_[cellI]))
        {
            bool interfaceCell = interface_[cellI];

            closestNormal = vector::zero;
            foundNormals.clear();

            loopStencilFindNearNeighbour
            (
                stencil[cellI],
                map,
                globalIdx,
                centre,
                normal,
                mesh_.C()[cellI],
                closestNormal,
                foundNormals
            );

            bool tooCoarse = false;

            if (foundNormals.size() > 1 && mag(closestNormal) != 0)
            {

                forAll(foundNormals, i)
                {
                    // all have the length of 1
                    // to coarse if normal angle is bigger than 45 deg
                    if ((closestNormal & foundNormals[i]) <= 0.707)
                    {
                        tooCoarse = true;
                    }
                }
            }
            else
            {
                tooCoarse = true;
            }

            // if a normal was found and the interface is fine enough
            // smallDist is always smallDist
            if (mag(closestNormal) != 0 && !tooCoarse)
            {
                surfnormals[cellI] = closestNormal;
            }
            else
            {
                surfnormals[cellI] = lsGrad.calcGrad(
                    cellI,mesh_.C()[cellI], alpha1_, stencil[cellI], mapAlpha, mapCC);
            }
        }
    }
}

void Foam::heightFunc::interpolateNormals
(
    const volVectorField& centre,
    const volVectorField& normal,
    const volVectorField& U,
    volVectorField& surfnormals
)
{
    scalar dt = mesh_.time().deltaTValue();

    leastSquareCell lsGrad(mesh_, nei_.globalNumbering());

    // create PtrList
    PtrList<volVectorField> interfaceData(2);
    // fill ptrlist
    interfaceData.set(0, centre);
    interfaceData.set(1, normal);

    // create maps
    Map<List<vector> > map;
    Map<scalar> mapAlpha;
    Map<vector> mapCC;

    // update maps and get stencil
    const List<DynamicList<label> >& stencil =
        nei_.getDataFromOtherProcs<scalar>(
            alpha1_, interface_, mapAlpha);
    // update maps
    nei_.getDataFromOtherProcs<vector>(mesh_.C(), interface_, mapCC);
    nei_.getDataFromOtherProcs<vector>(interfaceData, interface_, map);

    vector closestNormal = vector::zero;
    DynamicList<vector> foundNormals(10);

    const globalIndex& globalIdx = nei_.globalNumbering();

    forAll(alpha1_, cellI)
    {

        //if (isASurfaceCell(alpha1_[cellI]))
        if(interface_[cellI])
        {
            closestNormal = vector::zero;
            foundNormals.clear();

            loopStencilInterpolateNormal
            (
                stencil[cellI],
                map,
                globalIdx,
                centre,
                normal,
                mesh_.C()[cellI]-U[cellI]*dt,
                closestNormal,
                foundNormals
            );

            bool tooCoarse = false;

            if (foundNormals.size() > 1 && mag(closestNormal) != 0)
            {

                forAll(foundNormals, i)
                {
                    // all have the length of 1
                    // to coarse if normal angle is bigger than 10 deg
                    if ((closestNormal & foundNormals[i]) <= 0.98)
                    {
                        tooCoarse = true;
                    }
                }
            }
            else
            {
                tooCoarse = true;
            }

            // if a normal was found and the interface is fine enough
            // smallDist is always smallDist
            if (mag(closestNormal) != 0 && !tooCoarse)
            {
                surfnormals[cellI] = closestNormal;
            }
            else
            {
                surfnormals[cellI] = lsGrad.calcGrad(
                    cellI,mesh_.C()[cellI], alpha1_, stencil[cellI], mapAlpha, mapCC);
            }
        }
    }
}

// ************************************************************************* //
