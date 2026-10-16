Application
-----------

   phaseFieldFoam

Description
-----------

   Phase field solver based on work by Takada. This solver uses the
   Cahn-Hilliard equation and the Navier-Stokes coupling for the
   calculation of the phase field for two immiscible fluids by
   diffusive and advective transport mechanisms.

Original code by:
-----------

   Donaldson, Adam: Dalhousie University Halifax, Canada

Ported to OpenFOAM version 2.2.0 and further developed:
-------------------------------------------------------

   Weiss, Sebastian: TU Bergakademie Freiberg, Germany

Changes and fixes (2015-2018)
-----------------------------

Besides the port to OpenFOAM 2.2.0, the solver was corrected and extended
in several places. See the commit history for details.

**Physical corrections**

- Corrected the sign of the gravity/buoyancy contribution
  (`snGrad(rho)*ghf`, `phiHbyA += phig`).
- Corrected the sign of the surface force term (`surfaceForceTerm`).
- Removed the re-calculation of the face flux `phi` from the interpolated
  velocity, which overwrote the conservative flux.
- Mixture density and face viscosity are now evaluated consistently from the
  time-averaged phase field `tmpA = 0.5*(alpha1.oldTime() + alpha1)`
  (`rhoMix(tmpA)`, `muf(tmpA)`).

**Solver structure**

- Restructured the pressure-velocity coupling following the interFoam
  formalism (PIMPLE loop, `rhoPhi`, momentum equation, `UEqn.relax()`).
- Added a pre-conditioning step for the phase field and combined it with the
  solver (`preConditioner.H`, `preCAlphaEqn.H`, `alphaEqnSubCycle.H`).
- Added sub-cycling of the phase field equation.
- Added the turbulence contribution `turbulence->divDevRhoReff(rho, U)`
  to the momentum equation.
- Added adaptive time stepping (`setDeltaT.H`) based on the flow and
  interface Courant numbers.

**Dynamic meshes**

- Fixed the handling of moving meshes (`makeAbsolute(phi, U)`).
- Added a refinement indicator field `gradAlpha1Field` based on the
  capillary width and the phase field gradient, for adaptive mesh
  refinement at the interface.

**Robustness and output**

- Output of minimum phase field values on the boundaries.
- Replaced `double` with `scalar` for consistency with OpenFOAM's precision
  settings.
- `transportProperties` is re-read when modified during a run.

Further ports
-------------

   OpenFOAM 3.0.0: Weiss, Sebastian (phaseFieldFoam-300)
   OpenFOAM v6:    Segre, Enrico, Weizmann Institute, Israel (phaseFieldFoam-6)
