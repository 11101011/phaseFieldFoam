{
    scalar difCoNum = 
    (
        max
        (
            mesh.surfaceInterpolation::deltaCoeffs()*mag 
            (
                (scalar(1) + twoPhaseProperties.temperature())/Foam::pow(scalar(2),scalar(2)*twoPhaseProperties.temperature() + scalar(1))*gAlpha1
              + sqr(twoPhaseProperties.capillaryWidth())/scalar(4)*glapAlpha1
            )/mesh.magSf()
        ).value()*runTime.deltaT().value()*scalar(2)
    );

    difCoNum = max(difCoNum,scalar(1e-16));
    dimensionedScalar newMobility("mobility",dimensionSet(0,2,-1,0,0,0,0), twoPhaseProperties.mobilityCourantNo()/difCoNum);
    twoPhaseProperties.setMobility(newMobility);
}
