//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file homogeneous_sphere.cpp
//! \brief Problem generator for the homogeneous sphere problem.
//!
//! This problem calculates the radiative transfer through and outside of a sphere of 
//! constant density, temperature, and opacity.
//! Input parameters are:
//!    - problem/rho0    = density in the sphere
//!    - problem/T0      = temperature in the sphere 
//!    - problem/opac = opacity in the sphere (inverse length units)
//!
//! REFERENCE: Yan-Fei Jiang (姜燕飞), "An Implicit Finite Volume Scheme to Solve the 
//! Time-dependent Radiation Transport Equation Based on Discrete Ordinates",
//! ApJS, 253, 49 (2021), and references therein.
//========================================================================================

// C headers

// C++ headers
#include <iostream>   // endl
#include <fstream>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cmath>      // sqrt
#include <algorithm>  // min

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../hydro/hydro.hpp"
#include "../eos/eos.hpp"
#include "../bvals/bvals.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../field/field.hpp"
#include "../coordinates/coordinates.hpp"
#include "../nr_radiation/radiation.hpp"
#include "../nr_radiation/integrators/rad_integrators.hpp"
#include "../nr_radiation/implicit/radiation_implicit.hpp"

void RadConstantFluxInnerX1(
     MeshBlock *pmb, Coordinates *pco, NRRadiation *prad, 
     const AthenaArray<Real> &w, FaceField &b, 
     AthenaArray<Real> &ir,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void RadConstantFluxOuterX1(
     MeshBlock *pmb, Coordinates *pco, NRRadiation *prad, 
     const AthenaArray<Real> &w, FaceField &b, 
     AthenaArray<Real> &ir,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void HydroOuterX1(
    MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void HydroInnerX1(
    MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);



Real T0;
Real rho0;
Real opac;

void Mesh::InitUserMeshData(ParameterInput *pin) {

  T0 = pin->GetReal("problem", "T0");
  rho0 = pin->GetReal("problem", "rho0");
  opac = pin->GetReal("problem", "opac");

  // std::cout << "InitUserMeshData - (rho0,T0,kap): " << rho0 << " ; " << T0 << " ; "  << opac << std::endl;

  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, HydroInnerX1);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, HydroOuterX1);

  if(NR_RADIATION_ENABLED || IM_RADIATION_ENABLED) {
    // Enroll radiation boundary functions
    EnrollUserRadBoundaryFunction(BoundaryFace::inner_x1, RadConstantFluxInnerX1);
    EnrollUserRadBoundaryFunction(BoundaryFace::outer_x1, RadConstantFluxOuterX1);
  }

  return;
}


//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief homogeneous sphere test
//======================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
    
  Real er   = std::pow(T0, 4.0); 
  Real gamma = peos->GetGamma();

  // Initialize hydro variable
  for(int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        Real x1 = pcoord->x1v(i);
        if (x1 < 1.0){
            phydro->u(IDN,k,j,i) = rho0;
            phydro->u(IEN,k,j,i) = T0 * phydro->u(IDN,k,j,i)/(gamma-1.0);

        }
        else {
            phydro->u(IDN,k,j,i) = 1e-7;
            phydro->u(IEN,k,j,i) = 0.0;
        }
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
      }
    }
  }
  
  // Initialize opacity and specific intensity
  if(NR_RADIATION_ENABLED || IM_RADIATION_ENABLED){
    int nfreq = pnrrad->nfreq;

    Real *ir_lab;
    
    for(int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          
          ir_lab = &(pnrrad->ir(k,j,i,0));
          for(int n=0; n<pnrrad->n_fre_ang; n++){
             ir_lab[n] = er;
          }
          
          for (int ifr=0; ifr < nfreq; ++ifr){
            pnrrad->sigma_s(k,j,i,ifr) = 0.0;
            if (pcoord->x1v(i) < 1.0) {
                pnrrad->sigma_a(k,j,i,ifr) = opac;  // rosseland mean absorption opac
                pnrrad->sigma_pe(k,j,i,ifr) = opac; // radiation energy weighted mean
                pnrrad->sigma_p(k,j,i,ifr) = opac;  // planck mean
            }
            else {
                pnrrad->sigma_a(k,j,i,ifr) = 0.0;
                pnrrad->sigma_pe(k,j,i,ifr) = 0.0;
                pnrrad->sigma_p(k,j,i,ifr) = 0.0;
            }
          }
        }
      }
    }
    
  }// End Rad
  return;
}



//The following returns hydro variables to their initial values at each timestep
void MeshBlock::UserWorkInLoop(void)
{
    Real gamma = peos->GetGamma();
    for(int k=ks; k<=ke; k++) {
        for(int j=js; j<=je; j++) {
            for(int i=is; i<=ie; i++) {
                Real x1 = pcoord->x1v(i);
                if (x1 < 1.0){
                    phydro->u(IDN,k,j,i) = rho0;
                    phydro->u(IEN,k,j,i) = T0 * phydro->u(IDN,k,j,i)/(gamma-1.0);
                }
                else {
                    phydro->u(IDN,k,j,i) = 1e-7;
                    phydro->u(IEN,k,j,i) = 0.0;
                }
                phydro->u(IM1,k,j,i) = 0.0;
                phydro->u(IM2,k,j,i) = 0.0;
                phydro->u(IM3,k,j,i) = 0.0;
            }
        }
    }
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {

  // Not currently used
  if(NR_RADIATION_ENABLED) {
    // Enroll Opacity Function
    //prad->EnrollOpacityFunction(Opacity);
  }
  return;

}


// Inner Rad BC: Constant intensity isotropic
void RadConstantFluxInnerX1(
     MeshBlock *pmb, Coordinates *pco, NRRadiation *prad, 
     const AthenaArray<Real> &w, FaceField &b, 
     AthenaArray<Real> &ir,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh) {

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
        for(int ifr=0; ifr<prad->nfreq; ++ifr){
          for(int n=0; n<prad->nang; ++n){
            Real &miuz = prad->mu(0,k,j,is-i,n);
            Real &weight = prad->wmu(n);
            Real ris = pco->x1v(is);
            Real ris1 = pco->x1v(is+1);
            Real r= pco->x1v(is-i);
            Real didr = (prad->ir(k,j,is+1,ifr*prad->nang+n)-prad->ir(k,j,is,ifr*prad->nang+n))/
                         (ris1 - ris);

             // Set constant specific intensity of arT_0^4/4pi (=1 in dimensionless units)
             prad->ir(k,j,is-i,ifr*prad->nang+n) = 1.0;

          }
        }
      }
    }
  }
}

// Outer Rad BC: Constant intensity outward, zero inward
void RadConstantFluxOuterX1(
     MeshBlock *pmb, Coordinates *pco, NRRadiation *prad, 
     const AthenaArray<Real> &w, FaceField &b, 
     AthenaArray<Real> &ir,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh) {

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
        for(int ifr=0; ifr<prad->nfreq; ++ifr){
          for(int n=0; n<prad->nang; ++n){
            Real &miuz = prad->mu(0,k,j,ie+i,n);
            Real &weight = prad->wmu(n);
            if(miuz > 0.0){
              prad->ir(k,j,ie+i,ifr*prad->nang+n) = prad->ir(k,j,ie,ifr*prad->nang+n);
            } else {
              prad->ir(k,j,ie+i,ifr*prad->nang+n) = 0.;
            }
          }
        }
      }
    }
  }

}

// Inner BC: Constant density & pressure, zero velocity
void HydroInnerX1(
    MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{
  for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          prim(IDN,k,j,is-i) = rho0;
          prim(IVX,k,j,is-i) = 0.0;
          prim(IVY,k,j,is-i) = 0.0;
          prim(IVZ,k,j,is-i) = 0.0;
          prim(IPR,k,j,is-i) = rho0*T0;
        }
      }
    }
}

// Outer BC: very small density & pressure, zero velocity
void HydroOuterX1(
    MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          prim(IDN,k,j,ie+i) = 1.e-7;
          prim(IVX,k,j,ie+i) = 0.0;
          prim(IVY,k,j,ie+i) = 0.0;
          prim(IVZ,k,j,ie+i) = 0.0;
          prim(IPR,k,j,ie+i) = 1.e-12;
        }
      }
    }
}
