/*---------------------------------------------------------------------------*\
    Modified work | Copyright (c) 2017-2019, German Aerospace Center (DLR)
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


#include "isoAdvector.H"
#include "addToRunTimeSelectionTable.H"

#include "fvcSurfaceIntegrate.H"
#include "upwind.H"
#include "interpolationCellPoint.H"

#include "MULES.H"
#include "DynamicField.H"



// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace advection
{
    defineTypeNameAndDebug(isoAdvector, 0);
    addToRunTimeSelectionTable(advectionSchemes,isoAdvector, components);
}
}

template<typename Type>
Type Foam::advection::isoAdvector::faceValue
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& f,
    const label facei
) const
{
    if (mesh_.isInternalFace(facei))
    {
        return f.primitiveField()[facei];
    }
    else
    {
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

        // Boundary face. Find out which face of which patch
        const label patchi = pbm.patchID()[facei - mesh_.nInternalFaces()];

        if (patchi < 0 || patchi >= pbm.size())
        {
            FatalErrorInFunction
                << "Cannot find patch for face " << facei
                << abort(FatalError);
        }

        // Handle empty patches
        const polyPatch& pp = pbm[patchi];
        if (isA<emptyPolyPatch>(pp) || pp.empty())
        {
            return pTraits<Type>::zero;
        }

        const label patchFacei = pp.whichFace(facei);
        return f.boundaryField()[patchi][patchFacei];
    }
}


template<typename Type>
void Foam::advection::isoAdvector::setFaceValue
(
    GeometricField<Type, fvsPatchField, surfaceMesh>& f,
    const label facei,
    const Type& value
) const
{
    if (mesh_.isInternalFace(facei))
    {
        f.primitiveFieldRef()[facei] = value;
    }
    else
    {
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

        // Boundary face. Find out which face of which patch
        const label patchi = pbm.patchID()[facei - mesh_.nInternalFaces()];

        if (patchi < 0 || patchi >= pbm.size())
        {
            FatalErrorInFunction
                << "Cannot find patch for face " << facei
                << abort(FatalError);
        }

        // Handle empty patches
        const polyPatch& pp = pbm[patchi];
        if (isA<emptyPolyPatch>(pp) || pp.empty())
        {
            return;
        }

        const label patchFacei = pp.whichFace(facei);

        f.boundaryFieldRef()[patchi][patchFacei] = value;
    }
}

void Foam::advection::isoAdvector::timeIntegratedFlux()
{

    const scalar dt = mesh_.time().deltaT().value();

    // Create object for interpolating velocity to isoface centres
    interpolationCellPoint<vector> UInterp(U_);

    // For each downwind face of each surface cell we "isoadvect" to find dVf
    label nSurfaceCells = 0;

    // Clear out the data for re-use and reset list containing information
    // whether cells could possibly need bounding
    clearIsoFaceData();
    checkBounding_ = false;

    // Get necessary references
    const scalarField& phiIn = phi_.primitiveField();
    const scalarField& magSfIn = mesh_.magSf().primitiveField();


    // Get necessary mesh data
    const labelListList& CC = mesh_.cellCells();


    // Storage for isoFace points. Only used if writeIsoFacesToFile_
    DynamicList< List<point> > isoFacePts;
    const DynamicField<label>& interfaceLabels = surf_->interfaceLabels();

    // Loop through cells
    forAll(interfaceLabels, i)
    {
        const label celli = interfaceLabels[i];
        if (mag(surf_->normal()[celli]) != 0)
        {
             
            // This is a surface cell, increment counter, append and mark cell
            nSurfaceCells++;
            surfCells_.append(celli);


            // Cell is cut
                const point x0 = surf_->centre()[celli];
                vector n0 = -surf_->normal()[celli];
                n0 /= (mag(n0));

                // Get the speed of the isoface by interpolating velocity and
                // dotting it with isoface normal
                const scalar Un0 = UInterp.interpolate(x0, celli) & n0;
                const cell& cellFaces = mesh_.cells()[celli];

                forAll(cellFaces, fi)
                {
                    // Get current face index
                    const label facei = cellFaces[fi];

                    // Check if the face is internal face
                    if (mesh_.isInternalFace(facei))
                    {
                        bool isDownwindFace = false;
                        label otherCell = -1;

                        // Check if the cell is owner
                        if (celli == mesh_.owner()[facei])
                        {

                            if (phiIn[facei] > 10*SMALL)
                            {
                                isDownwindFace = true;
                            }

                            // Other cell is neighbour
                            otherCell = mesh_.neighbour()[facei];
                        }
                        else //celli is the neighbour
                        {
                            if (phiIn[facei] < -10*SMALL)
                            {
                                isDownwindFace = true;
                            }

                            // Other cell is the owner
                            otherCell = mesh_.owner()[facei];
                        }

                        // Calculate time integrated flux if face is downwind
                        if (isDownwindFace)
                        {
                            dVf_[facei] =  advectFace_.timeIntegratedFlux
                            (
                                   facei,
                                   x0,
                                   n0,
                                   Un0,
                                   dt,
                                   phiIn[facei],
                                   magSfIn[facei]
                            );
                        }

                        // We want to check bounding of neighbour cells to surface
                        // cells as well:
                        checkBounding_[otherCell] = true;

                        // Also check neighbours of neighbours.
                        // Note: consider making it a run time selectable
                        // extension level (easily done with recursion):
                        // 0 - only neighbours
                        // 1 - neighbours of neighbours
                        // 2 - ...
                        const labelList& nNeighbourCells = CC[otherCell];
                        forAll(nNeighbourCells, ni)
                        {
                            checkBounding_[nNeighbourCells[ni]] = true;
                        }
                    }
                    else
                    {
                        bsFaces_.append(facei);
                        bsx0_.append(x0);
                        bsn0_.append(n0);
                        bsUn0_.append(Un0);
                        // Note: we must not check if the face is on the
                        // processor patch here.
                    }
                }

        }
    }

    // Get references to boundary fields and mesh
    const surfaceScalarField::Boundary& phib = phi_.boundaryField();
    const surfaceScalarField::Boundary& magSfb = mesh_.magSf().boundaryField();
    surfaceScalarField::Boundary& dVfb = dVf_.boundaryFieldRef();

    //Is it really necessary to have both boundaryMesh and pBoundaryMesh?
    const fvBoundaryMesh& boundaryMesh = mesh_.boundary();
    const polyBoundaryMesh& pBoundaryMesh = mesh_.boundaryMesh();

    // Loop through boundary surface faces
    forAll(bsFaces_, fi)
    {
        // Get boundary face index (in the global list)
        const label facei = bsFaces_[fi];

        // Get necessary labels
        // Note: consider optimisation since whichPatch is expensive
        const label patchi = pBoundaryMesh.whichPatch(facei);
        const label start = boundaryMesh[patchi].patch().start();
        const label size = boundaryMesh[patchi].size();

        if (size > 0)
        {
            // Get patch local label
            const label patchFacei = facei - start;
            const scalar phiP = phib[patchi][patchFacei];

            if (phiP > 10*SMALL)
            {
                const scalar magSf = magSfb[patchi][patchFacei];

                dVfb[patchi][patchFacei] = advectFace_.timeIntegratedFlux
                (
                    facei,
                    bsx0_[fi],
                    bsn0_[fi],
                    bsUn0_[fi],
                    dt,
                    phiP,
                    magSf
                );



                // Check if the face is on processor patch and append it to
                // the list if necessary
                checkIfOnProcPatch(facei);
            }
        }
    }
	Info<< "Number of isoAdvector surface = "
	  << returnReduce(nSurfaceCells, sumOp<label>()) << endl;
}



void Foam::advection::isoAdvector::setDownwindFaces
(
    const label celli,
    DynamicLabelList& downwindFaces
) const
{
    DebugInFunction << endl;

    // Get necessary mesh data and cell information
    const labelList& own = mesh_.faceOwner();
    const cellList& cells = mesh_.cells();
    const cell& c = cells[celli];

    downwindFaces.clear();

    // Check all faces of the cell
    forAll(c, fi)
    {
        // Get face and corresponding flux
        const label facei = c[fi];
        const scalar phi = faceValue(phi_, facei);

        if (own[facei] == celli)
        {
            if (phi > 10*SMALL)
            {
                downwindFaces.append(facei);
            }
        }
        else if (phi < -10*SMALL)
        {
            downwindFaces.append(facei);
        }
    }

    downwindFaces.shrink();
}


void Foam::advection::isoAdvector::limitFluxes()
{

    // Get time step size
    const scalar dt = mesh_.time().deltaT().value();

//    scalarField alphaNew = alpha1In_ - fvc::surfaceIntegrate(dVf_);
    const scalar aTol = 1.0e-12;          // Note: tolerances
    const scalar maxAlphaMinus1 = 1;      // max(alphaNew - 1);
    const scalar minAlpha = -1;           // min(alphaNew);
    const label nUndershoots = 20;        // sum(neg0(alphaNew + aTol));
    const label nOvershoots = 20;         // sum(pos0(alphaNew - 1 - aTol));
    cellIsBounded_ = false;

    // Loop number of bounding steps
    for (label n = 0; n < nAlphaBounds_; n++)
    {
    Info<< "isoAdvector: bounding iteration " << n + 1 << endl;

        if (maxAlphaMinus1 > aTol) // Note: tolerances
        {
            surfaceScalarField dVfcorrected("dVfcorrected", dVf_);
            DynamicList<label> correctedFaces(3*nOvershoots);
            boundFromAbove(alpha1_.primitiveFieldRef(), dVfcorrected, correctedFaces);

            forAll(correctedFaces, fi)
            {
                label facei = correctedFaces[fi];

                // Change to treat boundaries consistently
                setFaceValue(dVf_, facei, faceValue(dVfcorrected, facei));
            }

            syncProcPatches(dVf_, phi_);
        }

        if (minAlpha < -aTol) // Note: tolerances
        {
            DebugInfo << "Bound from below... " << endl;

            scalarField alpha2 = 1.0 - alpha1_.primitiveFieldRef();
            surfaceScalarField dVfcorrected
            (
                "dVfcorrected",
                phi_*dimensionedScalar("dt", dimTime, dt) - dVf_
            );
            // If phi_ > 0 then dVf_ > 0 and mag(phi_*dt-dVf_) < mag(phi_*dt) as
            // it should.
            // If phi_ < 0 then dVf_ < 0 and mag(phi_*dt-dVf_) < mag(phi_*dt) as
            // it should.

            DynamicList<label> correctedFaces(3*nUndershoots);
            boundFromAbove(alpha2, dVfcorrected, correctedFaces);

            forAll(correctedFaces, fi)
            {
                label facei = correctedFaces[fi];

                // Change to treat boundaries consistently
                scalar phi = faceValue(phi_, facei);
                scalar dVcorr = faceValue(dVfcorrected, facei);
                setFaceValue(dVf_, facei, phi*dt - dVcorr);
            }

            syncProcPatches(dVf_, phi_);
        }
    }
}

void Foam::advection::isoAdvector::limitFluxes(const volScalarField::Internal& Sp,const volScalarField::Internal& Su)
{

    // Get time step size
    const scalar dt = mesh_.time().deltaT().value();

//    scalarField alphaNew = alpha1In_ - fvc::surfaceIntegrate(dVf_);
    const scalar aTol = 1.0e-12;          // Note: tolerances
    const scalar maxAlphaMinus1 = 1;      // max(alphaNew - 1);
    const scalar minAlpha = -1;           // min(alphaNew);
    const label nUndershoots = 20;        // sum(neg0(alphaNew + aTol));
    const label nOvershoots = 20;         // sum(pos0(alphaNew - 1 - aTol));
    cellIsBounded_ = false;

    // Loop number of bounding steps
    for (label n = 0; n < nAlphaBounds_; n++)
    {
	Info<< "isoAdvector: bounding iteration " << n + 1 << endl;
        if (maxAlphaMinus1 > aTol) // Note: tolerances always true!!!!
        {
            surfaceScalarField dVfcorrected("dVfcorrected", dVf_);
            DynamicList<label> correctedFaces(3*nOvershoots);
            boundFromAbove(alpha1_.ref(), dVfcorrected, correctedFaces,Sp,Su);

            forAll(correctedFaces, fi)
            {
                label facei = correctedFaces[fi];

                // Change to treat boundaries consistently
                setFaceValue(dVf_, facei, faceValue(dVfcorrected, facei));
            }

            syncProcPatches(dVf_, phi_);
        }

        if (minAlpha < -aTol) // Note: tolerances always true!!!!
        {
            DebugInfo << "Bound from below... " << endl;

           /* scalarField alpha2 = 1.0 - alpha1_.ref();
            surfaceScalarField dVfcorrected
            (
                "dVfcorrected",
                phi_*dimensionedScalar("dt", dimTime, dt) - dVf_
            );
            // If phi_ > 0 then dVf_ > 0 and mag(phi_*dt-dVf_) < mag(phi_*dt) as
            // it should.
            // If phi_ < 0 then dVf_ < 0 and mag(phi_*dt-dVf_) < mag(phi_*dt) as
            // it should.

            DynamicList<label> correctedFaces(3*nUndershoots);
            boundFromBelow(alpha2, dVfcorrected, correctedFaces,Sp,Su);

            forAll(correctedFaces, fi)
            {
                label facei = correctedFaces[fi];

                // Change to treat boundaries consistently
                scalar phi = faceValue(phi_, facei);
                scalar dVcorr = faceValue(dVfcorrected, facei);
                faceValue(dVf_, facei, phi*dt - dVcorr);
            }

            syncProcPatches(dVf_, phi_);*/

            surfaceScalarField dVfcorrected("dVfcorrected", dVf_);
            DynamicList<label> correctedFaces(3*nOvershoots);
            boundFromBelow(alpha1_.ref(), dVfcorrected, correctedFaces,Sp,Su);

            forAll(correctedFaces, fi)
            {
                label facei = correctedFaces[fi];

                // Change to treat boundaries consistently
                setFaceValue(dVf_, facei, faceValue(dVfcorrected, facei));
            }

            syncProcPatches(dVf_, phi_);
        }
    }
}


void Foam::advection::isoAdvector::boundFromAbove
(
    const scalarField& alpha1,
    surfaceScalarField& dVf,
    DynamicList<label>& correctedFaces
)
{

    // Get time step size
    const scalarField& meshV = mesh_.cellVolumes();
    const scalar dt = mesh_.time().deltaT().value();

    correctedFaces.clear();
    scalar aTol = 10*SMALL; // Note: tolerances
    scalar maxOvershoot = -GREAT;
    label maxOvershootCell = -1;

    DynamicList<label> downwindFaces(10);
    DynamicList<label> facesToPassFluidThrough(downwindFaces.size());
    DynamicList<scalar> dVfmax(downwindFaces.size());
    DynamicList<scalar> phi(downwindFaces.size());

    // Loop through alpha cell centred field
    forAll(alpha1, celli)
    {
        if (checkBounding_[celli])
        {
            const scalar Vi = meshV[celli];
            scalar alpha1New = alpha1[celli] - netFlux(dVf, celli)/Vi;
            scalar alphaOvershoot = alpha1New - 1.0;
            scalar fluidToPassOn = alphaOvershoot*Vi;
            label nFacesToPassFluidThrough = 1;

            if (alphaOvershoot > maxOvershoot)
            {
                maxOvershoot = alphaOvershoot;
                maxOvershootCell = celli;
            }

            bool firstLoop = true;
            // First try to pass surplus fluid on to neighbour cells that are
            // not filled and to which dVf < phi*dt
            while (alphaOvershoot > aTol && nFacesToPassFluidThrough > 0)
            {
                cellIsBounded_[celli] = true;

		facesToPassFluidThrough.clear();
	        dVfmax.clear();
                phi.clear();

                // Find potential neighbour cells to pass surplus phase to
                //DynamicList<label> downwindFaces(mesh_.cells()[celli].size());
                //getDownwindFaces(celli, downwindFaces);
		setDownwindFaces(celli, downwindFaces);

                //DynamicList<label> facesToPassFluidThrough(downwindFaces.size());
                //DynamicList<scalar> dVfmax(downwindFaces.size());
                //DynamicList<scalar> phi(downwindFaces.size());

                scalar dVftot = 0.0;
                nFacesToPassFluidThrough = 0;

                forAll(downwindFaces, fi)
                {
                    const label facei = downwindFaces[fi];
                    const scalar phif = faceValue(phi_, facei);
                    const scalar dVff = faceValue(dVf, facei);
                    const scalar maxExtraFaceFluidTrans = mag(phif*dt - dVff);

		    DebugInfo
                        << "downwindFace " << facei
                        << " has maxExtraFaceFluidTrans = "
                        << maxExtraFaceFluidTrans << endl;

                    if (maxExtraFaceFluidTrans/Vi > aTol)
//                    if (maxExtraFaceFluidTrans/Vi > aTol &&
//                    mag(dVfIn[facei])/Vi > aTol) //Last condition may be
//                    important because without this we will flux through uncut
//                    downwind faces
                    {
                        facesToPassFluidThrough.append(facei);
                        phi.append(phif);
                        dVfmax.append(maxExtraFaceFluidTrans);
                        dVftot += mag(phif*dt);
                    }
                }

                forAll(facesToPassFluidThrough, fi)
                {
                    const label facei = facesToPassFluidThrough[fi];
                    scalar fluidToPassThroughFace =
                        fluidToPassOn*mag(phi[fi]*dt)/dVftot;

                    nFacesToPassFluidThrough +=
                        pos0(dVfmax[fi] - fluidToPassThroughFace);

                    fluidToPassThroughFace =
                        min(fluidToPassThroughFace, dVfmax[fi]);

                    scalar dVff = faceValue(dVf, facei);
                    dVff += sign(phi[fi])*fluidToPassThroughFace;
                    setFaceValue(dVf, facei, dVff);

                    if (firstLoop)
                    {
                        checkIfOnProcPatch(facei);
                        correctedFaces.append(facei);
                    }
                }

                firstLoop = false;
                alpha1New = alpha1[celli] - netFlux(dVf, celli)/Vi;
                alphaOvershoot = alpha1New - 1.0;
                fluidToPassOn = alphaOvershoot*Vi;
		
		DebugInfo
                    << "\nNew alpha for cell " << celli << ": "
                    << alpha1New << endl;

            }
        }
    }
    DebugInfo << "correctedFaces = " << correctedFaces << endl;
}

void Foam::advection::isoAdvector::boundFromAbove
(
    const scalarField& alpha1,
    surfaceScalarField& dVf,
    DynamicList<label>& correctedFaces,
    const volScalarField::Internal& Sp,
    const volScalarField::Internal& Su
)
{

    // Get time step size
    const scalar dt = mesh_.time().deltaT().value();
    const scalar rdT = 1/mesh_.time().deltaT().value();

    correctedFaces.clear();
    scalar aTol = 10*SMALL; // Note: tolerances
    scalar maxOvershoot = -GREAT;
    label maxOvershootCell = -1;

    DynamicList<label> downwindFaces(10);
    DynamicList<label> facesToPassFluidThrough(downwindFaces.size());
    DynamicList<scalar> dVfmax(downwindFaces.size());
    DynamicList<scalar> phi(downwindFaces.size());

    // Loop through alpha cell centred field
    forAll(alpha1, celli)
    {
        if (checkBounding_[celli])
        {
            const scalar& Vi = mesh_.V()[celli];
            scalar alpha1New = (alpha1[celli]*rdT + Su[celli] - netFlux(dVf, celli)/Vi*rdT)/(rdT-Sp[celli]);
            scalar alphaOvershoot = alpha1New - 1.0;
            scalar fluidToPassOn = alphaOvershoot*Vi;
            label nFacesToPassFluidThrough = 1;

            if (alphaOvershoot > maxOvershoot)
            {
                maxOvershoot = alphaOvershoot;
                maxOvershootCell = celli;
            }

            bool firstLoop = true;
            // First try to pass surplus fluid on to neighbour cells that are
            // not filled and to which dVf < phi*dt
            while (alphaOvershoot > aTol && nFacesToPassFluidThrough > 0)
            {
                cellIsBounded_[celli] = true;

		facesToPassFluidThrough.clear();
	        dVfmax.clear();
                phi.clear();

                // Find potential neighbour cells to pass surplus phase to
                //DynamicList<label> downwindFaces(mesh_.cells()[celli].size());
                //getDownwindFaces(celli, downwindFaces);
		setDownwindFaces(celli, downwindFaces);

                //DynamicList<label> facesToPassFluidThrough(downwindFaces.size());
                //DynamicList<scalar> dVfmax(downwindFaces.size());
                //DynamicList<scalar> phi(downwindFaces.size());

                scalar dVftot = 0.0;
                nFacesToPassFluidThrough = 0;

                forAll(downwindFaces, fi)
                {
                    const label facei = downwindFaces[fi];
                    const scalar phif = faceValue(phi_, facei);
                    const scalar dVff = faceValue(dVf, facei);
                    const scalar maxExtraFaceFluidTrans = mag(phif*dt - dVff);

                    if (maxExtraFaceFluidTrans/Vi > aTol)
//                    if (maxExtraFaceFluidTrans/Vi > aTol &&
//                    mag(dVfIn[facei])/Vi > aTol) //Last condition may be
//                    important because without this we will flux through uncut
//                    downwind faces
                    {
                        facesToPassFluidThrough.append(facei);
                        phi.append(phif);
                        dVfmax.append(maxExtraFaceFluidTrans);
                        dVftot += mag(phif*dt);
                    }
                }

                forAll(facesToPassFluidThrough, fi)
                {
                    const label facei = facesToPassFluidThrough[fi];
                    scalar fluidToPassThroughFace =
                        fluidToPassOn*mag(phi[fi]*dt)/dVftot;

                    nFacesToPassFluidThrough +=
                        pos0(dVfmax[fi] - fluidToPassThroughFace);

                    fluidToPassThroughFace =
                        min(fluidToPassThroughFace, dVfmax[fi]);

                    scalar dVff = faceValue(dVf, facei);
                    dVff += sign(phi[fi])*fluidToPassThroughFace;
                    setFaceValue(dVf, facei, dVff);

                    if (firstLoop)
                    {
                        checkIfOnProcPatch(facei);
                        correctedFaces.append(facei);
                    }
                }

                firstLoop = false;
                alpha1New = (alpha1[celli]*rdT + Su[celli] - netFlux(dVf, celli)/Vi*rdT)/(rdT-Sp[celli]);
                alphaOvershoot = alpha1New - 1.0;
                fluidToPassOn = alphaOvershoot*Vi;

            }
        }
    }
}

void Foam::advection::isoAdvector::boundFromBelow
(
    const scalarField& alpha1,
    surfaceScalarField& dVf,
    DynamicList<label>& correctedFaces,
    const volScalarField::Internal& Sp,
    const volScalarField::Internal& Su
)
{

    // Get time step size
    const scalar dt = mesh_.time().deltaT().value();
    const scalar rdT = 1/mesh_.time().deltaT().value();

    correctedFaces.clear();
    scalar aTol = 10*SMALL; // Note: tolerances
    scalar maxOvershoot = -GREAT;
    label maxOvershootCell = -1;

    DynamicList<label> downwindFaces(10);
    DynamicList<label> facesToPassFluidThrough(downwindFaces.size());
    DynamicList<scalar> dVfmax(downwindFaces.size());
    //DynamicList<scalar> dV(downwindFaces.size());
    DynamicList<scalar> phi(downwindFaces.size());

    // Loop through alpha cell centred field
    forAll(alpha1, celli)
    {
        if (checkBounding_[celli])
        {
            const scalar& Vi = mesh_.V()[celli];
            scalar alpha1New = (alpha1[celli]*rdT + Su[celli] - netFlux(dVf, celli)/Vi*rdT)/(rdT-Sp[celli]);
            scalar alphaUndershoot = mag(alpha1New)*pos0(-alpha1New); // pos is with zero
            scalar fluidToHoldBack = alphaUndershoot*Vi;
            label nFacesToPassFluidThrough = 1;


            if (alphaUndershoot > maxOvershoot) // does make any sense is always true
            {
                maxOvershoot = alphaUndershoot; // why?
                maxOvershootCell = celli; // why?
            }

            bool firstLoop = true;
            // First try to pass surplus fluid on to neighbour cells that are
            // not filled and to which dVf < phi*dt

		    while (alphaUndershoot > aTol && nFacesToPassFluidThrough > 0)
		    {
		        cellIsBounded_[celli] = true;

			facesToPassFluidThrough.clear();
		        dVfmax.clear();
		        phi.clear();

		        // Find potential neighbour cells to pass surplus phase to
		        //DynamicList<label> downwindFaces(mesh_.cells()[celli].size());
		        //getDownwindFaces(celli, downwindFaces);
			setDownwindFaces(celli, downwindFaces);

		        //DynamicList<label> facesToPassFluidThrough(downwindFaces.size());
		        //DynamicList<scalar> dVfmax(downwindFaces.size());
		        //DynamicList<scalar> dV(downwindFaces.size());

		        scalar dVftot = 0.0;
		        nFacesToPassFluidThrough = 0;

		        forAll(downwindFaces, fi)
		        {
		            const label facei = downwindFaces[fi];
		            const scalar phif = faceValue(phi_, facei);
		            const scalar dVff = faceValue(dVf, facei);
		            //const scalar maxExtraFaceFluidTrans = mag(phif*dt - dVff);
			    const scalar maxFaceFluidTrans = mag(dVff);//mag(dVff);//mag( dVff); // zero is the lower limit phi the upper

		            if (maxFaceFluidTrans/Vi > aTol)
	//                    if (maxExtraFaceFluidTrans/Vi > aTol &&
	//                    mag(dVfIn[facei])/Vi > aTol) //Last condition may be
	//                    important because without this we will flux through uncut
	//                    downwind faces
		            {
		                facesToPassFluidThrough.append(facei);
                        phi.append(phif); // no it gets confuesing
		                dVfmax.append(maxFaceFluidTrans);
                        dVftot += mag(phif*dt); // cannot be zero cell wouldn't be below zero
		            }
		        }

		        forAll(facesToPassFluidThrough, fi)
		        {
		            const label facei = facesToPassFluidThrough[fi];
		            scalar fluidToPassThroughFace =
		                fluidToHoldBack*mag(phi[fi]*dt)/dVftot;

		            nFacesToPassFluidThrough +=
		                pos0(dVfmax[fi] - fluidToPassThroughFace);

		            fluidToPassThroughFace =
		                min(fluidToPassThroughFace, dVfmax[fi]);

		            scalar dVff = faceValue(dVf, facei);
		            dVff -= sign(phi[fi])*fluidToPassThroughFace;
		            setFaceValue(dVf, facei, dVff);

		            if (firstLoop)
		            {
		                checkIfOnProcPatch(facei);
		                correctedFaces.append(facei);
		            }
		        }

		        firstLoop = false;
		        alpha1New = (alpha1[celli]*rdT + Su[celli] - netFlux(dVf, celli)/Vi*rdT)/(rdT-Sp[celli]);
		        alphaUndershoot =  mag(alpha1New)*pos0(-alpha1New);
		
		        fluidToHoldBack = alphaUndershoot*Vi;

		    }
        }
    }
}

Foam::scalar Foam::advection::isoAdvector::netFlux
(
    const surfaceScalarField& dVf,
    const label celli
) const
{
    scalar dV = 0;

    // Get face indices
    const cell& c = mesh_.cells()[celli];

    // Get mesh data
    const labelList& own = mesh_.faceOwner();

    forAll(c, fi)
    {
        const label facei = c[fi];
        const scalar dVff = faceValue(dVf, facei);

        if (own[facei] == celli)
        {
            dV += dVff;
        }
        else
        {
            dV -= dVff;
        }
    }

    return dV;
}


void Foam::advection::isoAdvector::syncProcPatches
(
    surfaceScalarField& dVf,
    const surfaceScalarField& phi
)
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    if (Pstream::parRun())
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Send
        forAll(procPatchLabels_, i)
        {
            const label patchi = procPatchLabels_[i];

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchi]);

            UOPstream toNbr(procPatch.neighbProcNo(), pBufs);
            const scalarField& pFlux = dVf.boundaryField()[patchi];

            const List<label>& surfCellFacesOnProcPatch =
                surfaceCellFacesOnProcPatches_[patchi];

            const UIndirectList<scalar> dVfPatch
            (
                pFlux,
                surfCellFacesOnProcPatch
            );

            toNbr << surfCellFacesOnProcPatch << dVfPatch;
        }

        pBufs.finishedSends();


        // Receive and combine
        forAll(procPatchLabels_, patchLabeli)
        {
            const label patchi = procPatchLabels_[patchLabeli];

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchi]);

            UIPstream fromNeighb(procPatch.neighbProcNo(), pBufs);
            List<label> faceIDs;
            List<scalar> nbrdVfs;

            fromNeighb >> faceIDs >> nbrdVfs;

            if (debug)
            {
                Pout<< "Received at time = " << mesh_.time().value()
                    << ": surfCellFacesOnProcPatch = " << faceIDs << nl
                    << "Received at time = " << mesh_.time().value()
                    << ": dVfPatch = " << nbrdVfs << endl;
            }

            // Combine fluxes
            scalarField& localFlux = dVf.boundaryFieldRef()[patchi];

            forAll(faceIDs, i)
            {
                const label facei = faceIDs[i];
                localFlux[facei] = - nbrdVfs[i];
                if (debug && mag(localFlux[facei] + nbrdVfs[i]) > 10*SMALL)
                {
                    Pout<< "localFlux[facei] = " << localFlux[facei]
                        << " and nbrdVfs[i] = " << nbrdVfs[i]
                        << " for facei = " << facei << endl;
                }
            }
        }

        if (debug)
        {
            // Write out results for checking
            forAll(procPatchLabels_, patchLabeli)
            {
                const label patchi = procPatchLabels_[patchLabeli];
                const scalarField& localFlux = dVf.boundaryField()[patchi];
                Pout<< "time = " << mesh_.time().value() << ": localFlux = "
                    << localFlux << endl;
            }
        }

        // Reinitialising list used for minimal parallel communication
        forAll(surfaceCellFacesOnProcPatches_, patchi)
        {
            surfaceCellFacesOnProcPatches_[patchi].clear();
        }
    }
}


void Foam::advection::isoAdvector::checkIfOnProcPatch(const label facei)
{
    if (!mesh_.isInternalFace(facei))
    {
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
        const label patchi = pbm.patchID()[facei - mesh_.nInternalFaces()];

        if (isA<processorPolyPatch>(pbm[patchi]) && pbm[patchi].size())
        {
            const label patchFacei = pbm[patchi].whichFace(facei);
            surfaceCellFacesOnProcPatches_[patchi].append(patchFacei);
        }
    }
}



void Foam::advection::isoAdvector::calcFaceFlux()
{
    const scalar dt = mesh_.time().deltaT().value();
    // Create object for interpolating velocity to isoface centres
    interpolationCellPoint<vector> UInterp(U_);

    const labelList& owner = mesh_.faceOwner();
    const labelList& neighbour = mesh_.faceNeighbour();

    forAll(dVf_,faceI) // all are internal
    {
        const label own = owner[faceI];
        const label nei = neighbour[faceI];
        {

            scalar weight = pos0(phi_[faceI]);
            if(weight == 1 && surf_().interfaceCell()[own])
            {
                // take own value
                vector n0 = surf_().normal()[own];
                n0/= mag(n0);

                // Get the speed of the isoface by interpolating velocity and
                // dotting it with isoface normal

                const scalar Un0 = UInterp.interpolate(surf_().centre()[own], own) & n0;

               dVf_[faceI] =  advectFace_.timeIntegratedFlux
               (
                   faceI,
                   surf_().centre()[own],
                   n0,
                   Un0,
                   dt,
                   phi_[faceI],
                   mesh_.magSf()[faceI]
               );
            }
            else if(weight == 0 &&  surf_().interfaceCell()[nei])
            {
                // take nei values
                vector n0 = surf_().normal()[nei];
                n0/= mag(n0);

                // Get the speed of the isoface by interpolating velocity and
                // dotting it with isoface normal

                const scalar Un0 = UInterp.interpolate(surf_().centre()[nei], nei) & n0;

                dVf_[faceI] =  advectFace_.timeIntegratedFlux
                (
                       faceI,
                       surf_().centre()[nei],
                       n0,
                       Un0,
                       dt,
                       phi_[faceI],
                       mesh_.magSf()[faceI]
                );
            }
            else
            {
                dVf_[faceI] = phi_[faceI]*(weight*(alpha1_[own] - alpha1_[nei]) + alpha1_[nei])*dt;
            }
        }
    }

    surfaceScalarField::Boundary& dVfb = dVf_.boundaryFieldRef();
    const surfaceScalarField::Boundary& phib = phi_.boundaryField();

    forAll(dVfb,patchI) // loop patches
    {
        forAll(dVfb[patchI],i) // for now proc patches
        {
            scalar weight = pos0(phib[patchI][i]);
            if(weight == 1) // downwind face
            {

            }
            else
            {
                dVfb[patchI][i] = alpha1_.boundaryField()[patchI][i]*phib[patchI][i]*dt;
            }
         }
    }
}

void Foam::advection::isoAdvector::applyBruteForceBounding()
{
    bool alpha1Changed = false;

    scalar snapAlphaTol = modelDict().lookupOrDefault<scalar>("snapTol", 0);
    if (snapAlphaTol > 0)
    {
        alpha1_ =
            alpha1_
           *pos0(alpha1_ - snapAlphaTol)
           *neg0(alpha1_ - (1.0 - snapAlphaTol))
          + pos0(alpha1_ - (1.0 - snapAlphaTol));

        alpha1Changed = true;
    }

    bool clip = modelDict().lookupOrDefault<bool>("clip", true);
    if (clip)
    {
        alpha1_ = min(scalar(1.0), max(scalar(0.0), alpha1_));
        alpha1Changed = true;
    }

    if (alpha1Changed)
    {
        alpha1_.correctBoundaryConditions();
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::advection::isoAdvector::isoAdvector
(
        volScalarField& alpha1,
        const surfaceScalarField& phi,
        const volVectorField& U
)
:
    advectionSchemes
    (
        typeName,
        alpha1,
        phi,
        U
    ),
    mesh_(alpha1.mesh()),
    dVf_
    (
        IOobject
        (
            "dVf_",
            alpha1.mesh().time().timeName(),
            alpha1.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha1.mesh(),
        dimensionedScalar("zero", dimVol, 0)
    ),
    advectFace_(alpha1.mesh(),alpha1),
    // Tolerances and solution controls
    nAlphaBounds_(modelDict().lookupOrDefault<label>("nAlphaBounds", 3)),
    vof2IsoTol_(modelDict().lookupOrDefault<scalar>("vof2IsoTol", 1e-8)),
    surfCellTol_(modelDict().lookupOrDefault<scalar>("surfCellTol", 1e-8)),
    gradAlphaBasedNormal_(modelDict().lookupOrDefault<bool>("gradAlphaNormal", false)),
    writeIsoFacesToFile_(modelDict().lookupOrDefault<bool>("isoFaces2File", false)),

    // Cell cutting data
    surfCells_(label(0.2*mesh_.nCells())),
    cellIsBounded_(mesh_.nCells(), false),
    checkBounding_(mesh_.nCells(), false),
    bsFaces_(label(0.2*(mesh_.nFaces() - mesh_.nInternalFaces()))),
    bsx0_(bsFaces_.size()),
    bsn0_(bsFaces_.size()),
    bsUn0_(bsFaces_.size()),
    minMagSf_(gMin(mesh_.magSf())),
    limiter_(0),

    // Parallel run data
    procPatchLabels_(mesh_.boundary().size()),
    surfaceCellFacesOnProcPatches_(0)
{

    // Prepare lists used in parallel runs
    if (Pstream::parRun())
    {
        // Force calculation of required demand driven data (else parallel
        // communication may crash)
        mesh_.cellCentres();
        mesh_.cellVolumes();
        mesh_.faceCentres();
        mesh_.faceAreas();
        mesh_.magSf();
        mesh_.boundaryMesh().patchID();
        mesh_.cellPoints();
        mesh_.cellCells();
        mesh_.cells();

        // Get boundary mesh and resize the list for parallel comms
        const polyBoundaryMesh& patches = mesh_.boundaryMesh();

        surfaceCellFacesOnProcPatches_.resize(patches.size());

        // Append all processor patch labels to the list
        forAll(patches, patchi)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchi])
             && patches[patchi].size() > 0
            )
            {
                procPatchLabels_.append(patchi);
            }
        }
    }

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::advection::isoAdvector::~isoAdvector()
{}

// * * * * * * * * * * * * * * Protected Access Member Functions  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //
void Foam::advection::isoAdvector::advect()
{
    mesh_.time().cpuTimeIncrement();
    // set to zero counts the time between calls
    // only on master?
    surf_->reconstruct();
    runTime_.x() += mesh_.time().cpuTimeIncrement();


    // Initialising dVf with upwind values
    dVf_ = upwind<scalar>(mesh_, phi_).flux(alpha1_)*mesh_.time().deltaT();

    // Do the isoAdvector on surface cells
    // Reconstruct advect
    timeIntegratedFlux();

    // Synchronize processor patches
    syncProcPatches(dVf_, phi_);

    // Adjust dVf for unbounded cells
    limitFluxes();


    // Advect the free surface
    alpha1_ -= fvc::surfaceIntegrate(dVf_);
    alpha1_.correctBoundaryConditions();

    applyBruteForceBounding();

    alphaPhi_ = dVf_/mesh_.time().deltaT();

    runTime_.y() += mesh_.time().cpuTimeIncrement();

    reduce(runTime_.x(),maxOp<scalar>());
    reduce(runTime_.y(),maxOp<scalar>());
    Info << "runTime " << runTime_ << endl;

}

void Foam::advection::isoAdvector::advect(const volScalarField::Internal& Sp,const volScalarField::Internal& Su)
{
// splitted in div part and source term part 
// source term part get calculated first
    scalar rDeltaT = 1/mesh_.time().deltaTValue();

// MULES similiar to MULES
//    alpha1_.primitiveFieldRef() =
//    (
//      alpha1_.primitiveField()*rDeltaT
//      + Su.field()
//      - fvc::surfaceIntegrate(dVf_)
//    )/(rDeltaT - Sp.field());


// second part
    surf_->reconstruct();

//    calcFaceFlux();
    // Initialising dVf with upwind values
    dVf_ = upwind<scalar>(mesh_, phi_).flux(alpha1_)*mesh_.time().deltaT();

    // Do the isoAdvector on surface cells
    // Reconstruct advect
    timeIntegratedFlux();

    // Synchronize processor patches
    syncProcPatches(dVf_, phi_);

   /* alpha1_.primitiveFieldRef() =
    (
      alpha1_.primitiveField()*rDeltaT
      + Su.field()
      - fvc::surfaceIntegrate(dVf_)().primitiveField()*rDeltaT
    )/(rDeltaT - Sp.field());
    alpha1_.correctBoundaryConditions();*/

    // Adjust dVf for unbounded cells
    limitFluxes(Sp,Su);
// MULES

//   psi.primitiveFieldRef() =
//   (
//      rho.field()*psi.primitiveField()*rDeltaT
//      + Su.field()
//      - psiIf
//   )/(rho.field()*rDeltaT - Sp.field());

    alpha1_.primitiveFieldRef() =
    (
      alpha1_.oldTime().primitiveField()*rDeltaT
      + Su.field()
      - fvc::surfaceIntegrate(dVf_)().primitiveField()*rDeltaT
    )/(rDeltaT - Sp.field());

  //  alpha1_ -= fvc::surfaceIntegrate(dVf_) ; // dvf mal deltaT
    alpha1_.correctBoundaryConditions();

    applyBruteForceBounding();

    alphaPhi_ = dVf_/mesh_.time().deltaT();

}

