/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    twoPhaseMixture

Description
    This class contains functions for calculating the following fields based on
    the current alpha1 field: density, viscosity, mobility, derivative of 
    Mixing Energy
    The class is also responsible for calculating constants specific to the 
    phase field model, such as kappaThickness, kappaSurfaceT

\*---------------------------------------------------------------------------*/

#include "twoPhaseMixture.H"
#include "addToRunTimeSelectionTable.H"
#include "dimensionedScalar.H"
#include "fixedGradientFvPatchFields.H"
#include "surfaceFields.H"
#include "fvc.H"
#include "mathematicalConstants.H"


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //
namespace Foam
{
using namespace constant::mathematical;

scalar twoPhaseMixture::calc2F1(double a1_,double a2_, double b_, double z_)
{
    scalar result = 1;
    scalar c_ = 1;
    for (int iCounter=0; iCounter<10000; iCounter++)
    {
        c_ *= (iCounter + a1_)*(iCounter + a2_)/(iCounter + scalar(1))/(iCounter + b_)*z_;
        result += c_;
    }
    return result; 
}

//-Calculate and return the laminar viscosity
void twoPhaseMixture::calcNu()
{
    nuModel1_->correct();
    nuModel2_->correct();

    const volScalarField limitedAlpha1
    (
        "limitedAlpha1",
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    //-Average kinematic viscosity calculated from dynamic viscosity
    nu_ = mu(alpha1_)/(limitedAlpha1*rho1_ + (scalar(1) - limitedAlpha1)*rho2_);
}

//-Calculate the capillary width which controls interfacial thickness
void twoPhaseMixture::calcCapillaryWidth()
{
    //-This formulation applies to the TVSED energy function
    //-Variable units: [capillaryWidth] = (0 1 0 0 0 0 0)

    capillaryWidth_ = thickness_/Foam::pow(scalar(2),Tr_.value() + scalar(0.5))/(scalar(1) 
  - scalar(2)*filterAlpha_)/calc2F1(scalar(0.5),(scalar(1) + Tr_.value())/scalar(2),scalar(1.5),sqr(scalar(1)
  - scalar(2)*filterAlpha_.value()));

    Info<< "Capillary Width = " 
        << capillaryWidth_.value() 
        << endl;
}

//-Calculates the mixing energy density, which controls the surface tension
void twoPhaseMixture::calcMixingEDensity()
{
    //-The following entry is for the standard surface tension expression
    mixingEDensity_ = 
    (
        sigma_*capillaryWidth_*Foam::pow(scalar(2.0),scalar(0.5)
      + Tr_.value())/calc2F1(scalar(0.5),(Tr_.value()
      + scalar(1))/scalar(-2),scalar(1.5),scalar(1))
    );

    Info<< "Mixing Energy Density = " 
        << mixingEDensity_.value() 
        << endl;
}

//-Returns the TVSED Energy function value at curAlpha1_
dimensionedScalar twoPhaseMixture::mixingE(const scalar curAlpha1_)
{
    dimensionedScalar scaledAlpha1_ ("alpha1",dimensionSet(0,0,0,0,0,0,0),curAlpha1_);
    return Foam::pow(mag(scaledAlpha1_*(scalar(1) - scaledAlpha1_)),scalar(1) + Tr_.value());
}

//-Slope function for the contact angle gradient calculation
scalarField twoPhaseMixture::boundarySlope(const scalarField& curAlpha1_)
{
    scalarField slope_ =
    (
        scalar(-1)*cos
        (
            theta_.value()*pi/scalar(180)
        )*sqrt
        (
            scalar(2)*mixingEscalar(curAlpha1_)
        )/capillaryWidth_.value()
    );

    return slope_;
}

scalarField twoPhaseMixture::mixingEscalar(const scalarField& cAlpha1_)
{
    scalarField limitedAlpha1 = min(max(cAlpha1_,scalar(0)),scalar(1));

    return tmp<scalarField> 
    (
        new scalarField
        (
            Foam::pow(mag(limitedAlpha1*(scalar(1) - limitedAlpha1)), scalar(1) + Tr_.value())
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

twoPhaseMixture::twoPhaseMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const IOdictionary& dict,
    const word& alpha1Name
)
:
    transportModel(U, phi),

    phase1Name_(found("phases") ? wordList(lookup("phases"))[0] : "phase1"),
    phase2Name_(found("phases") ? wordList(lookup("phases"))[1] : "phase2"),

    nuModel1_
    (
        viscosityModel::New
        (
            "nu1",
            subDict(phase1Name_),
            U,
            phi
        )
    ),
    nuModel2_
    (
        viscosityModel::New
        (
            "nu2",
            subDict(phase2Name_),
            U,
            phi
        )
    ),

    rho1_("rho", dimDensity, nuModel1_->viscosityProperties().lookup("rho")),
    rho2_("rho", dimDensity, nuModel2_->viscosityProperties().lookup("rho")),

    filterAlpha_ (dict.lookup("filterAlpha")),
    thickness_ (dict.lookup("thickness")),
    Tr_ (dict.lookup("Tr")),
    sigma_(dict.lookup("sigma")),
    mobilityCourant_(dict.lookup("mobilityCourant")),
    theta_(dict.lookup("theta")),

    U_(U),
    phi_(phi),
 
    alpha1_
    (
        IOobject
        (
            found("phases") ? word("alpha" + phase1Name_) : alpha1Name,
            U_.time().timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh()
    ),

    nu_
    (
        IOobject
        (
            "nu",
            U_.time().timeName(),
            U_.db()
        ),
        U_.mesh(),
        dimensionedScalar("nu", dimensionSet(0,2,-1,0,0),0),
        calculatedFvPatchScalarField::typeName
    ),

    mobility_
    (
        "mobility",
        dimensionSet(0,2,-1,0,0,0,0),
        scalar(0)
    ),

    mixingEDensity_
    (
        "mixingEDensity",
        dimensionSet(1,1,-2,0,0,0,0), 
        scalar(0)
    ),

    capillaryWidth_
    (
        "capillaryWidth",
        dimensionSet(0,1,0,0,0,0,0),
        scalar(0)
    )

{
    calcNu();
    calcCapillaryWidth();
    calcMixingEDensity();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
tmp<volScalarField> twoPhaseMixture::mu(const volScalarField& alpha1New_) const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1New_,scalar(0)),scalar(1))
    );

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "mu",
            limitedAlpha1*rho1_*nuModel1_->nu()
          + (scalar(1) - limitedAlpha1)*rho2_*nuModel2_->nu()
        )
    );
}


tmp<surfaceScalarField> twoPhaseMixture::muf(const volScalarField& alpha1New_) const
{
    const surfaceScalarField alpha1f
    (
        min(max(fvc::interpolate(alpha1New_),scalar(0)),scalar(1))
    );

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "muf",
            alpha1f*rho1_*fvc::interpolate(nuModel1_->nu()) + (scalar(1) - alpha1f)*rho2_*fvc::interpolate(nuModel2_->nu())
        )
    );
}


tmp<Foam::surfaceScalarField> Foam::twoPhaseMixture::nuf() const
{
    const surfaceScalarField alpha1f
    (
        min(max(fvc::interpolate(alpha1_),scalar(0)),scalar(1))
    );

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "nuf",
            (
                alpha1f*rho1_*fvc::interpolate(nuModel1_->nu()) 
              + (scalar(1) - alpha1f)*rho2_*fvc::interpolate(nuModel2_->nu())
            )/(alpha1f*rho1_ + (scalar(1) - alpha1f)*rho2_)
        )
    );
}


tmp<surfaceScalarField> twoPhaseMixture::diffusivityF(const surfaceScalarField& alpha1New_) const
{
    scalar safetyFactor = scalar(0);
    if (Tr_.value() < scalar(0.001))
    {
        safetyFactor = scalar(0.001);
    }

    return tmp<surfaceScalarField> 
    (
        new surfaceScalarField
        (
            "diffusivity",
            (scalar(1) + Tr_.value())*
            (
                Tr_.value()
              - (
                    scalar(4)*Tr_.value() + scalar(2)
                )*alpha1New_
            )*Foam::pow
            (
                mag
                (
                    alpha1New_
                )
              + safetyFactor
              , Tr_.value()
            )
        )
    );
}

tmp<volScalarField> twoPhaseMixture::rhoMix(const volScalarField& alpha1New_) const
{
    volScalarField limitedAlpha1 = (min(max(alpha1New_,scalar(0)),scalar(1)));
    return tmp<volScalarField> 
    (
        new volScalarField
        (
            "rho",
            (scalar(1) - limitedAlpha1)*rho2_ + limitedAlpha1*rho1_
        )
    );
}

tmp<surfaceScalarField> twoPhaseMixture::rhoMixF(const volScalarField& alpha1New_) const
{
    surfaceScalarField limitedAlpha1f = (min(max(fvc::interpolate(alpha1New_),scalar(0)),scalar(1)));

    return tmp<surfaceScalarField> 
    (
        new surfaceScalarField
        (
            "rho",
            (
                (scalar(1) - limitedAlpha1f)*rho2_ + limitedAlpha1f*rho1_
            )
        )
    );
}

//-return angle term for alpha1 to density conversion
tmp<volScalarField> twoPhaseMixture::alpha1Angle(const volScalarField& alpha1New_) const
{
    return tmp<volScalarField> 
    (
        new volScalarField
        (
            "angle",
            (min(max(alpha1New_,filterAlpha_),scalar(1) - filterAlpha_) - scalar(0.5))/(scalar(1) - scalar(2)*filterAlpha_)*pi
        )
    );
}

//-Return multiplier term
dimensionedScalar twoPhaseMixture::alpha1Multiplier() const
{
    return sqr(scalar(0.5)*pi/(scalar(1) - scalar(2)*filterAlpha_));
}

void twoPhaseMixture::updateContactAngle(volScalarField& curAlpha1_)
{
    //-Create a pointer to the mesh
    const fvMesh& mesh = curAlpha1_.mesh();

    //-Obtain a list of all boundaries on the mesh
    const fvPatchList& patches = mesh.boundary();

    //-Cycle through each boundary, current boundary indicated by patchi within the loop
    forAll(patches, patchi)
    {
        //-Check to see if the current boundary is a fixedGradient type
        if (isA<fixedGradientFvPatchScalarField>(curAlpha1_.boundaryField()[patchi]))
        {
            //-Create a reference to the patch field. Note that this variable can be treated
            // as a scalar field containing the values of curAlpha1_ on the boundary face
            fixedGradientFvPatchScalarField& curPatch = refCast<fixedGradientFvPatchScalarField>(curAlpha1_.boundaryField()[patchi]);

            //-Create a pointer reference to the gradient field
            scalarField& gradAlpha1 = curPatch.gradient();

            //-Set the gradient to zero to get the adjacent cell values during the next evaluation call
            gradAlpha1 *= scalar(0);

            //-Ensure the value of alpha1 is updated on the boundary face.  Equivalent to correctBoundaryConditions function
            curPatch.evaluate();

            scalarField zeroAlpha1 = curPatch;

            //-Determine the new gradient value. The function "boundarySlope" accepts a scalarField, and returns a scalarField
            gradAlpha1 = boundarySlope(curPatch);

            //-Recalculate the value of curAlpha1 on the boundary faces
            curPatch.evaluate();

            scalarField upperFilter = pos(curPatch - scalar(1));
            scalarField lowerFilter = pos(scalar(-1)*curPatch);

            //-Create a safetyFactor to avoid division by zero when the boundary alpha1 = 0 or 1, or when a zero-gradient is used
            // Determines location where alpha1 is 0 & 1
            scalarField safetyFactor = (pos(curPatch - zeroAlpha1 + scalar(0.00001)) - neg(zeroAlpha1 - curPatch + scalar(0.00001)))*scalar(0.0001);

            //-Perform the slope adjustment to ensure the boundary value of alpha1 is >= 0 and <= 1.
            gradAlpha1 *= 
            (
                upperFilter*(scalar(1) - zeroAlpha1)/(curPatch - zeroAlpha1 + safetyFactor)
              + lowerFilter*(scalar(-1)*zeroAlpha1)/(curPatch - zeroAlpha1 + safetyFactor)
              + (scalar(1) - upperFilter - lowerFilter)
            );

            curPatch.evaluate();  

            Info<< "Boundary Min: " 
                << min(curPatch*(scalar(1) - curPatch)) 
                << endl;
        }
    }

    //-Simply ensure that all boundary conditions are updated (i.e. if you have other boundaries, such as zeroGradient, this is required)
    curAlpha1_.correctBoundaryConditions();
}

dimensionedScalar twoPhaseMixture::epsTOne() 
{
    //-This formulation applies to the TVSED energy function
    // Variable units: [capillaryWidth] = (0 1 0 0 0 0 0)
    return  thickness_/Foam::pow(scalar(2),scalar(1.5))/(scalar(1) - scalar(2)*filterAlpha_)/calc2F1(scalar(0.5),scalar(1),scalar(1.5),sqr(scalar(1)
  - scalar(2)*filterAlpha_.value()));
}

dimensionedScalar twoPhaseMixture::mixingEDensityTOne()
{
    return 
    (
        sigma_*Foam::pow(epsTOne(),scalar(2))*Foam::pow(scalar(2),scalar(2.5)
      - Tr_.value())/capillaryWidth_/(calc2F1(scalar(0.5),(Tr_.value()
      - scalar(3))/scalar(2),scalar(1.5),scalar(1)))
    );
}

bool twoPhaseMixture::read()
{
    if (transportModel::read())
    {
        if
        (
            nuModel1_().read(subDict(phase1Name_)) 
         && nuModel2_().read(subDict(phase2Name_))
        )
        {
            nuModel1_->viscosityProperties().lookup("rho") >> rho1_;
            nuModel2_->viscosityProperties().lookup("rho") >> rho2_;

            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}//-End namespace Foam

// ************************************************************************* //
