    //-Create a temporal alpha1
    volScalarField tmpA = (alpha1.oldTime() + alpha1)/scalar(2);

    //-Pre-calculate grad(alpha1) to reduce comp. time
    gradAlpha1 = fvc::grad(min(max(tmpA,scalar(0)),scalar(1)));

    //-Surface force using Kim's model
    volScalarField ICurv =
    fvc::div
    (
        (
            fvc::interpolate
            (
                gradAlpha1/(mag(gradAlpha1)
              + (scalar(1E-8)/pow(average(tmpA.mesh().V()),scalar(1)/scalar(3))))
            )
        )& mesh.Sf()
    );

    volVectorField surfaceForceTerm = - twoPhaseProperties.mixingEDensity()*ICurv*mag(gradAlpha1)*gradAlpha1;
    surfaceScalarField muf = twoPhaseProperties.muf(tmpA);

    //-Calculate the final uncorrected velocity field
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho,U)
      + fvm::div(rhoPhi,U)
      - fvm::laplacian(muf,U)
      - (fvc::grad(U)& fvc::grad(muf))
      + turbulence->divDevRhoReff(rho,U)
    );

    UEqn.relax();

    surfaceScalarField KSaG =
    (
      - fvc::snGrad(rho)*ghf*mesh.magSf()
      + (fvc::interpolate(surfaceForceTerm)& mesh.Sf())
      * neg(fvc::interpolate(tmpA*(tmpA - scalar(1))))
    );