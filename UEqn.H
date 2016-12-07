    //-Create a temporal alpha1
    tmpA = (alpha1 + alpha1.oldTime())/scalar(2);

    //-Pre-calculate grad(alpha1) to reduce comp. time
    gradAlpha1 = fvc::grad(min(max(tmpA,scalar(0)),scalar(1)));

    //-Surface force using Kim's model
    ICurv = fvc::div
    (
        (
            fvc::interpolate
            (
                gradAlpha1/(mag(gradAlpha1) + (scalar(1E-8)/pow(average(tmpA.mesh().V()),scalar(1)/scalar(3))))
            )
        )& mesh.Sf()
    );

    surfaceForceTerm = twoPhaseProperties.mixingEDensity()*ICurv*mag(gradAlpha1)*gradAlpha1;

    //-Calculate the final uncorrected velocity field
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho,U)
      + fvm::div(rhoPhi,U)
      - fvm::laplacian(twoPhaseProperties.muf(tmpA),U)
      - (fvc::grad(U)& fvc::grad(twoPhaseProperties.muf(tmpA)))
    );

    UEqn.relax();

    KSaG = 
    (
        fvc::snGrad(rho)*ghf*mesh.magSf()
      + (fvc::interpolate(surfaceForceTerm)& mesh.Sf())
      * neg(fvc::interpolate(tmpA*(tmpA - scalar(1))))
    );