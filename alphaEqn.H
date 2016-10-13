    {
        //-Solve the alpha1 equation at the current time step and
        // calculate the laplacian of the alpha1 field and the diffusive flux.

        surfaceScalarField gAlpha1 = fvc::snGrad(alpha1)*mesh.magSf();
        surfaceScalarField glapAlpha1 = fvc::snGrad(fvc::laplacian(alpha1))*mesh.magSf();

        #include "mobEqn.H"

        surfaceScalarField alpha1Energy = fvc::interpolate(alpha1)*(scalar(1) - fvc::interpolate(alpha1));

        difFlux = twoPhaseProperties.mobility()*pos(alpha1Energy)*
        (   
            (
                alpha1Energy*sqr(twoPhaseProperties.capillaryWidth())*glapAlpha1
            )
          - (
                twoPhaseProperties.diffusivityF(alpha1Energy)*gAlpha1
            )
        );

        #include "limitFlux.H"
    }

    volScalarField tempK_Alpha1
    (
        "tempK_Alpha1",
        (
            - (fvc::div(difFlux))
            - (fvc::div(phi, alpha1))
        )
    );

    rhoPhi = (difFlux + fvc::flux(phi,alpha1,"div(phi,alpha1)"))*(twoPhaseProperties.rho1() - twoPhaseProperties.rho2()) + phi*twoPhaseProperties.rho2();
    rhoPhiSum += rhoPhi*K_Multiplier;

    //-Solve for the new value of alpha1
    alpha1 = alpha1.oldTime() + runTime.deltaT()*T_Multiplier*tempK_Alpha1;

    twoPhaseProperties.updateContactAngle(alpha1);

    rho = twoPhaseProperties.rhoMix(scalar(0.5)*(alpha1 + alpha1.oldTime()));
    rhoPhi = rhoPhiSum;