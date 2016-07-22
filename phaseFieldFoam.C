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

Application
    phaseFieldFoam

Description
    Phase field solver based on work by Takada. This solver uses the 
    Cahn-Hilliard equation and the Navier-Stokes coupling for the 
    calculation of the phase field for two immiscible fluids by 
    diffusive and advective transport mechanisms.
    
    Written by:
    Donaldson, Adam: Dalhousie University Halifax, Canada
    
    Ported to OpenFOAM version 2.2.0:
    Weiss, Sebastian: TU Bergakademie Freiberg, Germany

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "interfaceProperties.H"
#include "twoPhaseMixture.H"
#include "pimpleControl.H"
#include "fvIOoptionList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    
    pimpleControl pimple(mesh);
    
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "readTimeControls.H"
    #include "correctPhi.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;
        
        //---------------------------------------------------//
        
        scalar timeBeforeMeshUpdate = runTime.elapsedCpuTime();
 
        // Do any mesh changes
        mesh.update();
    
        if (mesh.changing())
        {
            Info<< "Execution time for mesh.update() = "
                << runTime.elapsedCpuTime() - timeBeforeMeshUpdate
                << " s" << endl;

            gh = g& mesh.C();
            ghf = g& mesh.Cf();
        }

        if (mesh.changing() && correctPhi)
        {
            #include "correctPhi.H"
        }

        if (mesh.changing() && checkMeshCourantNo)
        {
            #include "meshCourantNo.H"
        }
        
        //---------------------------------------------------//

        fvc::makeRelative(phi, U);
        twoPhaseProperties.correct();
        
        //- Update the refinement field indicator
        gradAlpha1Field = twoPhaseProperties.capillaryWidth()*mag(fvc::grad(alpha1))
                        / Foam::pow(scalar(2),scalar(0.5))/Foam::pow(twoPhaseProperties.filterAlpha()*(scalar(1)
                        - twoPhaseProperties.filterAlpha()),(scalar(1) + twoPhaseProperties.temperature())*scalar(0.5));

        {
            volScalarField checkAlpha1 = scalar(10)*(pos(alpha1 - twoPhaseProperties.filterAlpha()/scalar(2)) - neg(scalar(1) 
                                       - twoPhaseProperties.filterAlpha()/scalar(2) - alpha1));

            gradAlpha1Field += checkAlpha1;
        }

        //- RungeKutta 4th order method
        volScalarField K_alpha1 ("K_alpha1", alpha1*scalar(0)/runTime.deltaT());
        surfaceScalarField rhoPhiSum = scalar(0)*rhoPhi;

        scalar T_Multiplier = scalar(0);
        scalar K_Multiplier = scalar(0);

        Info<< "Solving U and Alpha1 RK4 Equations: ";

        for (int i=0; i<=3; i++)
        {
            Info << " " << scalar(i);
            T_Multiplier = scalar(0.5) + scalar(i/2)*scalar(0.5);
            K_Multiplier = scalar(1)/(scalar(3) + scalar(3)*mag(scalar(1) - scalar((i + 1)/scalar(2))));

            #include "alphaEqn.H"
            K_alpha1 += K_Multiplier*tempK_Alpha1;
        }
        
        Info<< " ... Complete" << endl;

        alpha1 = alpha1.oldTime() + runTime.deltaT()*K_alpha1;
        twoPhaseProperties.updateContactAngle(alpha1);

        {
           volScalarField tempVolFrac = pos(alpha1 - scalar(0.5));
        
            Info<< "Phase-1 volume fraction = "
                << alpha1.weightedAverage(mesh.Vsc()).value()
                << "  Min(alpha1) = " << min(alpha1).value()
                << "  Max(alpha1) = " << max(alpha1).value()
                << endl;
        }
//         rho = twoPhaseProperties.rhoMix(scalar(0.5)*(alpha1 + alpha1.oldTime()));
        rho = alpha1*rho1 + (scalar(1) - alpha1)*rho2;
        rhoPhi = rhoPhiSum;

        //- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"

            //- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info << "End\n" << endl;
}

// ************************************************************************* //