if (pimple.nCorrPIMPLE() > 1)
{
    // If nOuterCorrectors > 1 then for all but the first loop the advection
    // of alpha is done using an average, 0.5*phi+0.5*phiNew where phi is
    // the flux at the beginning of the time step and phiNew is the flux
    // estimate at the end of the time step from the previous outer
    // iteration. Similarly we use 0.5*U + 0.5*UNew in later iterations.
    if (pimple.firstIter())
    {
        // To recalculate the alpha1 update in subsequent iterations, we
        // must store its current value before overwriting with the new
        // value
        alpha1.storePrevIter();
        // Storing initial phi and U for use in later outer iterations.
        phi.storePrevIter();
        U.storePrevIter();
    }
    else
    {
        // Resetting alpha1 to value before advection in first PIMPLE
        // iteration.
        alpha1 = alpha1.prevIter();

        // Setting U and phi with which to advect interface.
        U = 0.5*U.prevIter() + 0.5*U;
        phi = 0.5*phi.prevIter() + 0.5*phi;
    }
}

tmp<surfaceScalarField> talphaPhi1(alphaPhi10);

if (nAlphaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();

    talphaPhi1 = new surfaceScalarField
    (
        IOobject
        (
            "alphaPhi1",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(alphaPhi10.dimensions(), Zero)
    );

    surfaceScalarField rhoPhiSum
    (
        IOobject
        (
            "rhoPhiSum",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(rhoPhi.dimensions(), Zero)
    );

    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
        #include "alphaEqn.H"
        talphaPhi1.ref() += (runTime.deltaT()/totalDeltaT)*alphaPhi10;
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
    }

    rhoPhi = rhoPhiSum;
}
else
{
    #include "alphaEqn.H"
}

if (!pimple.firstIter())
{
    // Resetting U and phi to value at latest iteration.
    U = 2.0*U - U.prevIter();
    phi = 2.0*phi - phi.prevIter();
}

rho == alpha1*rho1 + alpha2*rho2;


const surfaceScalarField& alphaPhi1 = talphaPhi1();
surfaceScalarField alphaPhi2("alphaPhi2", phi - alphaPhi1);

volScalarField::Internal contErr
(
    (
        fvc::ddt(rho) + fvc::div(rhoPhi)
      - (fvOptions(alpha1, mixture.thermo1().rho())&rho1)
      - (fvOptions(alpha2, mixture.thermo2().rho())&rho2)
    )()
);
