{
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

// Updating alpha1
advector->advect();
#include "rhofs.H"
rhoPhi = advector->getRhoPhi(rho1f, rho2f);

    if (!pimple.firstIter())
    {
        // Resetting U and phi to value at latest iteration.
        U = 2.0*U - U.prevIter();
        phi = 2.0*phi - phi.prevIter();
    }

    alpha2 = 1.0 - alpha1;
    mixture.correct();

}

Info<< "Phase-1 volume fraction = "
    << alpha1.weightedAverage(mesh.Vsc()).value()
    << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
    << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
    << endl;
