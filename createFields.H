    Info<< "Reading field fraction of liquid phase\n" << endl;
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    surfaceScalarField difFlux
    (
        IOobject
        (
            "difFlux",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    #include "createPhi.H"
    
    //-For input, to avoid precalculation, the alpha1 field simply must
    // be = 1 for liquid, and 0 for gas, and vary in between to set the interface

    twoPhaseMixture twoPhaseProperties(U, phi, transportProperties);
    volScalarField& alpha1(twoPhaseProperties.alpha1());

    const dimensionedScalar& rho1 = twoPhaseProperties.rho1();
    const dimensionedScalar& rho2 = twoPhaseProperties.rho2();

    volVectorField gradAlpha1
    (
        IOobject
        (
            "grad(alpha1)",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::grad(alpha1)
    );

    //-Construct interface from alpha1 distribution
    interfaceProperties interface(alpha1, U, twoPhaseProperties);
    
    #include "readGravitationalAcceleration.H"

    phi = fvc::interpolate(U)& mesh.Sf();

    //-Need to store rho for ddt(rho, U)
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        twoPhaseProperties.rhoMix(alpha1),
        alpha1.boundaryField().types()
    );
    rho.oldTime();

    Info << "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g& mesh.C());
    surfaceScalarField ghf("ghf",g& mesh.Cf());
    
    volScalarField gradAlpha1Field
    (
        IOobject
        (
            "gradAlpha1Field",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        twoPhaseProperties.capillaryWidth()*mag(gradAlpha1)/Foam::pow(scalar(2),scalar(0.5))/Foam::pow(twoPhaseProperties.filterAlpha()*(scalar(1)
      - twoPhaseProperties.filterAlpha()),(scalar(1) + twoPhaseProperties.temperature())*scalar(0.5))
    );

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rho*gh
    );
    
    surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rho*phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(rho)*phi
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        p_rgh,
        mesh.solutionDict().subDict("PIMPLE"),
        pRefCell,
        pRefValue
    );

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }

    fv::IOoptionList fvOptions(mesh);
