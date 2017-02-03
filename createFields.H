
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

    Info<< "Reading/calculating face flux field phi\n" << endl;

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(U)& mesh.Sf()
    );


    //-For input, to avoid precalculation, the alpha1 field simply must
    // be = 1 for liquid, and 0 for gas, and vary in between to set the interface

    Info<< "Reading transportProperties\n" << endl;
    twoPhaseMixture twoPhaseProperties(U,phi,transportProperties);

    volScalarField& alpha1(twoPhaseProperties.alpha1());

    volScalarField rho = twoPhaseProperties.rhoMix(alpha1);
    rho.oldTime();

    scalar T_Multiplier = scalar(0);
    scalar K_Multiplier = scalar(0);

    const scalar deltaTZero = 1E-8;


    volScalarField tempAlpha1
    (
        IOobject
        (
            "tempAlpha1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha1.oldTime() - fvc::div(difFlux)*runTime.deltaT()*T_Multiplier
    );

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

    volScalarField nDelta
    (
        IOobject
        (
            "nDelta",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
      - fvc::div(fvc::interpolate(lowFilter)*difFlux)*runTime.deltaT()*T_Multiplier
    );

    surfaceScalarField gAlpha1
    (
        IOobject
        (
            "gAlpha1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::snGrad(alpha1)*mesh.magSf()
    );

    volScalarField tempK_Alpha1
    (
        IOobject
        (
            "tempK_Alpha1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
      - fvc::div(difFlux) - fvc::div(phi,alpha1)
    );

    surfaceScalarField glapAlpha1
    (
        IOobject
        (
            "glapAlpha1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::snGrad(fvc::laplacian(alpha1))*mesh.magSf()
    );

    surfaceScalarField alpha1Energy
    (
        IOobject
        (
            "alpha1Energy",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(alpha1)*(scalar(1) - fvc::interpolate(alpha1))
    );

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

    volScalarField K_alpha1
    (
        IOobject 
        (
            "K_alpha1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha1*scalar(0)/runTime.deltaT()
    );

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
      - twoPhaseProperties.filterAlpha()),(scalar(1)
      + twoPhaseProperties.temperature())*scalar(0.5))
    );

    //-Mass flux
    // Initialisation does not matter because rhoPhi is reset after the
    // alpha1 solution before it is used in the U equation.
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

    surfaceScalarField rhoPhiSum
    (
        IOobject
        (
            "rhoPhiSum",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        scalar(0)*rhoPhi
    );

    //-Construct interface from alpha1 distribution
    interfaceProperties interface(alpha1,U,twoPhaseProperties);

    //-Construct incompressible turbulence model
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, twoPhaseProperties)
    );

    #include "readGravitationalAcceleration.H"

    Info << "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g& mesh.C());
    surfaceScalarField ghf("ghf",g& mesh.Cf());

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
