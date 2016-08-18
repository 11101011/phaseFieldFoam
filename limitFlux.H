{
    volScalarField tempAlpha1 = alpha1.oldTime() - fvc::div(difFlux)*runTime.deltaT()*T_Multiplier;

    //-Limit alpha field at 0
    volScalarField lowFilter
    (
        IOobject
        (
            "lowFilter",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        neg(tempAlpha1)
    );

    volScalarField newDelta = - fvc::div(fvc::interpolate(lowFilter)*difFlux)*runTime.deltaT()*T_Multiplier;
    lowFilter *= min(max((alpha1.oldTime() + newDelta)*lowFilter/sign(newDelta)/max(mag(newDelta),scalar(1e-32)),scalar(-4)),scalar(5));

    //-Limit alpha field at 1
    volScalarField highFilter
    (
        IOobject
        (
            "highFilter",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        neg(scalar(1) - tempAlpha1)
    );

    newDelta = - fvc::div(fvc::interpolate(highFilter)*difFlux)*runTime.deltaT()*T_Multiplier;
    highFilter *= min(max((scalar(1) - alpha1.oldTime() - newDelta)*scalar(-1)*highFilter/sign(newDelta)/max(mag(newDelta),scalar(1e-32)),scalar(-4)),scalar(5));

    //- Update difFlux
    difFlux *= (scalar(1) - fvc::interpolate(highFilter) - fvc::interpolate(lowFilter));
}
