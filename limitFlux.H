{
    tempAlpha1 = alpha1.oldTime() - fvc::div(difFlux)*runTime.deltaT()*T_Multiplier;

    lowFilter = neg(tempAlpha1);

    nDelta = - fvc::div(fvc::interpolate(lowFilter)*difFlux)*runTime.deltaT()*T_Multiplier;

    lowFilter *= min(max((alpha1.oldTime() + nDelta)*lowFilter/sign(nDelta)/max(mag(nDelta),scalar(1E-32)),scalar(-4)),scalar(5));

    highFilter = neg(scalar(1) - tempAlpha1);

    nDelta = - fvc::div(fvc::interpolate(highFilter)*difFlux)*runTime.deltaT()*T_Multiplier;

    highFilter *= min(max((scalar(1) - alpha1.oldTime() - nDelta)*scalar(-1)*highFilter/sign(nDelta)/max(mag(nDelta),scalar(1E-32)),scalar(-4)),scalar(5));

    //-Update difFlux
    difFlux *= (scalar(1) - fvc::interpolate(highFilter) - fvc::interpolate(lowFilter));
}