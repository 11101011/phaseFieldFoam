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

SourceFiles
    twoPhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhaseMixture_H
#define twoPhaseMixture_H

#include "incompressible/transportModel/transportModel.H"
#include "incompressible/viscosityModels/viscosityModel/viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class twoPhaseMixture Declaration                      
\*---------------------------------------------------------------------------*/

class twoPhaseMixture
:
    public transportModel
{
protected:
    
    //-Private data

        word phase1Name_;
        word phase2Name_;

        autoPtr<viscosityModel> nuModel1_;
        autoPtr<viscosityModel> nuModel2_;

        dimensionedScalar rho1_;
        dimensionedScalar rho2_;

        dimensionedScalar filterAlpha_;
        dimensionedScalar thickness_;
        dimensionedScalar Tr_;
        dimensionedScalar sigma_;
        dimensionedScalar mobilityCourant_;
        dimensionedScalar theta_;

        const volVectorField& U_;
        const surfaceScalarField& phi_;

        volScalarField alpha1_;

        volScalarField nu_;

        dimensionedScalar mobility_;

        dimensionedScalar mixingEDensity_;
        dimensionedScalar capillaryWidth_;

    //-Private Member Functions

        scalar calc2F1(double a1_, double a2_, double b_, double z_);

        //-Calculate and return the laminar viscosity
        void calcNu();

        //-Calculate and record the capillary width
        void calcCapillaryWidth();

        //-Calculate and record the mixing Energy Density * Beta
        void calcMixingEDensity();

        //-Calculate and return the TVSED value at a given curAlpha1
        dimensionedScalar mixingE(const scalar curAlpha1_);

        //-Correct the alpha1 for a given contact angle
        scalarField boundarySlope(const scalarField& curAlpha1_);

        //-Calculate the mixing energy
        scalarField mixingEscalar(const scalarField& curAlpha1_);


public:

    //-Constructors

        //-Construct from components
        twoPhaseMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const IOdictionary& dict,
            const word& alpha1Name = "alpha1"
        );


    //-Destructor

    ~twoPhaseMixture()
    {}


    //-Member Functions

        const word phase1Name() const
        {
            return phase1Name_;
        }

        const word phase2Name() const
        {
            return phase2Name_;
        }

        scalar temperature()
        {
            return Tr_.value();
        }

        //-Return the phase-fraction of phase 1
        const volScalarField& alpha1() const
        {
            return alpha1_;
        }

        //-Return the phase-fraction of phase 1
        volScalarField& alpha1()
        {
            return alpha1_;
        }

        //-Return const-access to liquid phase viscosityModel
        const viscosityModel& nuModel1() const
        {
            return nuModel1_();
        }

        //-Return const-access to gas phase viscosityModel
        const viscosityModel& nuModel2() const
        {
            return nuModel2_();
        }

        //-Return const-access to the liquid density
        const dimensionedScalar& rho1() const
        {
            return rho1_;
        }

        //-Return const-access to gas density
        const dimensionedScalar& rho2() const
        {
            return rho2_;
        };

        //-Return const-access to liquid limit alpha1
        const dimensionedScalar& filterAlpha() const
        {
            return filterAlpha_;
        }

        //-Return the current mobility courant number
        const scalar& mobilityCourantNo() const
        {
            return mobilityCourant_.value();
        }

        const dimensionedScalar& mobility() const
        {
            return mobility_;
        }

        //-Set the mobilty 
        void setMobility(const dimensionedScalar newMob_)
        {
            mobility_ = newMob_;
        }

        //-Return the dynamic laminar viscosity
        tmp<volScalarField> mu(const volScalarField& alpha1New_) const;

        //-Return the face-interpolated dynamic laminar viscosity
        tmp<surfaceScalarField> muf(const volScalarField& alpha1New_) const;

        //-Return the kinematic laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //-Return the face-interpolated dynamic laminar viscosity
        tmp<surfaceScalarField> nuf() const;

        //-Return the diffusivity term
        tmp<surfaceScalarField> diffusivityF(const surfaceScalarField& alpha1New_) const;

        //-Return the kappaSurface
        dimensionedScalar mixingEDensity() const
        {
            return mixingEDensity_;
        }

        //-Return the kappaThickness
        dimensionedScalar capillaryWidth() const
        {
            return capillaryWidth_;
        }

        //-Return the mixture density
        tmp<volScalarField> rhoMix(const volScalarField& alpha1New_) const;

        //-Return the mixture density
        tmp<surfaceScalarField> rhoMixF(const volScalarField& alpha1New_) const;

        //-Return the alpha1 angle
        tmp<volScalarField> alpha1Angle(const volScalarField& alpha1New_) const;

        //-Return multiplier term
        dimensionedScalar alpha1Multiplier() const;

        //-Correct the laminar viscosity 
        virtual void correct()
        {
            calcNu();
        }

        //-Read base transportProperties dictionary
        virtual bool read();

        //-Correct the alpha1 for a given contact angle
        void updateContactAngle(volScalarField& curAlpha1_);

        dimensionedScalar epsTOne();
        dimensionedScalar mixingEDensityTOne();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}//-End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
