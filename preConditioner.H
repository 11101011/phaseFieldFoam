
    //-Update the refinement field indicator
    gradAlpha1Field =
    (
        twoPhaseProperties.capillaryWidth()*mag(fvc::grad(alpha1))
      / Foam::pow(scalar(2),scalar(0.5))/Foam::pow(twoPhaseProperties.filterAlpha()*(scalar(1)
      - twoPhaseProperties.filterAlpha()),(scalar(1) 
      + twoPhaseProperties.temperature())*scalar(0.5))
    );

    {
        volScalarField checkAlpha1 = 
        (
            scalar(10)*(pos(alpha1 
          - twoPhaseProperties.filterAlpha()/scalar(2)) 
          - neg(scalar(1) 
          - twoPhaseProperties.filterAlpha()/scalar(2) 
          - alpha1))
        );

        gradAlpha1Field += checkAlpha1;
    }

    //-Estimate relative flux
    fvc::makeRelative(phi,U);
    twoPhaseProperties.correct();

    //-RungeKutta 4th order method
    K_alpha1 = alpha1*scalar(0)/runTime.deltaT();

    rhoPhiSum = scalar(0)*rhoPhi;

    T_Multiplier = scalar(0);
    K_Multiplier = scalar(0);

    Info<< "Solving for alpha1 w/ Runge-Kutta algorithm: ";
    for (int i=0; i<=3; i++)
    {
        Info << " " << scalar(i);
        T_Multiplier = scalar(0.5) + scalar(i/2)*scalar(0.5);
        K_Multiplier = scalar(1)/(scalar(3) + scalar(3)*mag(scalar(1) - scalar((i + 1)/scalar(2))));

        #include "preCAlphaEqn.H"

        K_alpha1 += K_Multiplier*tempK_Alpha1;
    }

    Info << " ... Complete" << endl;

    alpha1 = alpha1.oldTime() + runTime.deltaT()*K_alpha1;

    twoPhaseProperties.updateContactAngle(alpha1);

    rho = twoPhaseProperties.rhoMix(scalar(0.5)*(alpha1.oldTime() + alpha1));

    rhoPhi = rhoPhiSum;

    #include "continuityErrs.H"
