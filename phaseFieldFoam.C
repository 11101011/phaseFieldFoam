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
#include "turbulenceModel.H"
#include "pimpleControl.H"
#include "fvIOoptionList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "readTimeControls.H"

    pimpleControl pimple(mesh);

    #include "correctPhi.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    bool t = true;
    bool b = true;

    //-Obtain a list of all boundaries on the mesh
    const fvPatchList& patches = mesh.boundary();
    int N = static_cast<int> (patches.size());

    scalar boundaryMin[N];
    bool boundaryMin_t[N];

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readControls.H"

        if (t && (runTime.value()<=runTime.deltaTValue()))
        {
            #include "CourantNo.H"
            #include "alphaCourantNo.H"

            //-Adjust time step for preCPhaseFieldFoam
            runTime.setDeltaT(deltaTZero);
        }
        else
        {
            #include "CourantNo.H"
            #include "alphaCourantNo.H"
            
            fvc::makeAbsolute(phi,U);
            
            #include "setDeltaT.H"
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        //---------------------------------------------------//

        //-Common part of phaseFieldFoam and preCPhaseFieldFoam
        scalar timeBeforeMeshUpdate = runTime.elapsedCpuTime();

        {
            //-Calculate the relative velocity used to map the relative flux phi
            volVectorField Urel("Urel", U);

            if (mesh.moving())
            {
                Urel -= fvc::reconstruct(fvc::meshPhi(U));
            }

            //-Update the refinement field indicator
            gradAlpha1Field = 
            twoPhaseProperties.capillaryWidth()*mag(fvc::grad(alpha1))/Foam::pow(scalar(2),scalar(0.5))/Foam::pow(twoPhaseProperties.filterAlpha()*(scalar(1)
          - twoPhaseProperties.filterAlpha()),(scalar(1)
          + twoPhaseProperties.temperature())*scalar(0.5));
                
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

            //-Do any mesh changes
            mesh.update();
        }

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

        if (t && (runTime.value()<=runTime.deltaTValue()))
        {
            Info<< nl << "Running pre-conditioner:" << nl << endl;
            t = false;

            #include "preConditioner.H"
        }
        else
        {
            //-After the first time step run the pressure & U loop
            if (b)
            {
                Info<< nl << "Running phase field calculation:" << nl << endl;
                b = false;
            }

            #include "alphaEqnSubCycle.H"

            //-Pressure-velocity PIMPLE corrector loop
            while (pimple.loop())
            {
                #include "UEqn.H"

                //- Pressure corrector loop
                while (pimple.correct())
                {
                    #include "pEqn.H"
                }
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;
}

// ************************************************************************* //
