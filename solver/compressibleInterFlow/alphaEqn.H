// If there are more than one outer corrector, we use a mixture of old and
// new U and phi for propagating alpha1 in all but the first outer iteration
if (!pimple.firstIter())
{
    // We are recalculating alpha1 from the its old time value
    alpha1 = alpha1.oldTime();
    // Temporarily storing new U and phi values in prevIter storage
    U.storePrevIter();
    phi.storePrevIter();

    // Overwriting new U and phi values with mixture of old and new values
    phi = 0.5*(phi + phi.oldTime());
    U = 0.5*(U + U.oldTime());
}

// Update alpha1
#include "alphaSuSp.H"
advector->advect(Sp,(Su + divU*min(alpha1(), scalar(1)))());
//advector->advect();

// Update rhoPhi
rhoPhi = advector->getRhoPhi(rho1, rho2);
alphaPhi10 = advector->alphaPhi();

alpha2 = 1.0 - alpha1;

if (!pimple.firstIter())
{
    // Restoring new U and phi values temporarily saved in prevIter() above
    U = U.prevIter();
    phi = phi.prevIter();
}

Info<< "Phase-1 volume fraction = "
    << alpha1.weightedAverage(mesh.Vsc()).value()
    << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
    << "  Max(" << alpha1.name() << ") - 1 = " << max(alpha1).value() - 1
    << endl;
