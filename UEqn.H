    //- Create a temporal alpha1
    volScalarField tmpA("alpha1", (alpha1 + alpha1.oldTime())/scalar(2));

    //- Pre-calculate grad(alpha1) to reduce comp. time
    gradAlpha1 = fvc::grad(min(max(tmpA,scalar(0)),scalar(1)));

    //- Surface force using kim's model
    volScalarField Curvature =
    fvc::div
    (
        (
            fvc::interpolate
            (
                gradAlpha1/(mag(gradAlpha1) + scalar(1e-8)/pow(average(tmpA.mesh().V()),scalar(1)/scalar(3)))
            )
        )& mesh.Sf()
    );

    volVectorField surfaceForceTerm
    (
        "div(lamdastar*beta*grad(alpha1)*grad(alpha1))",
        twoPhaseProperties.mixingEDensity()*Curvature*mag(gradAlpha1)*gradAlpha1
    );

    //-In the annex, muf is a function of tmpA
    surfaceScalarField muf = twoPhaseProperties.muf(tmpA);

    //- Calculate the final uncorrected velocity field
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      - fvm::laplacian(muf,U)
      - (fvc::grad(U)& fvc::grad(muf))
      + fvm::div(rhoPhi,U)
    );

    UEqn.relax();

    surfaceScalarField K_SandG
    (
        "K_SandG",
        (
          - fvc::snGrad(rho)*ghf*mesh.magSf()
          + (fvc::interpolate(surfaceForceTerm)& mesh.Sf())*neg(fvc::interpolate(tmpA*(tmpA - scalar(1))))
        )
    );

