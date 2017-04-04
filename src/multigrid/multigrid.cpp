//========================================================================================
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================

// C/C++ headers
#include <iostream>
#include <cmath>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cstring>    // memset

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "./multigrid.hpp"


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::LoadFinestData(const AthenaArray<Real> &src, int ns, int ngh)
//  \brief Fill the active zone of the finest level
void Multigrid::LoadFinestData(const AthenaArray<Real> &src, int ns, int ngh)
{
  AthenaArray<Real> &dst=u_[nlev_-1];
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+nx_-1, je=js+ny_-1, ke=ks+nz_-1;
#pragma ivdep
  for(int n=0; n<nvar; n++) {
    int nsrc=ns+n;
#pragma ivdep
    for(int k=ngh, mk=ks; mk<=ke; k++, mk++) {
#pragma ivdep
      for(int j=ngh, mj=js; mj<=je; j++, mj++) {
#pragma ivdep
        for(int i=ngh, mi=is; mi<=ie; i++, mi++)
          dst(n,mk,mj,mi)=src(nsrc,k,j,i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::LoadSource(const AthenaArray<Real> &src, int ns, int ngh,
//                                  Real fac)
//  \brief Fill the active zone of the finest level
void Multigrid::LoadSource(const AthenaArray<Real> &src, int ns, int ngh, Real fac)
{
  AthenaArray<Real> &dst=src_[nlev_-1];
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+nx_-1, je=js+ny_-1, ke=ks+nz_-1;
#pragma ivdep
  for(int n=0; n<nvar; n++) {
    int nsrc=ns+n;
#pragma ivdep
    for(int k=ngh, mk=ks; mk<=ke; k++, mk++) {
#pragma ivdep
      for(int j=ngh, mj=js; mj<=je; j++, mj++) {
#pragma ivdep
        for(int i=ngh, mi=is; mi<=ie; i++, mi++)
          dst(n,mk,mj,mi)=src(nsrc,k,j,i)*fac;
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::RetrieveResult(AthenaArray<Real> &dst, int ns, int ngh)
//  \brief Set the result, including the ghost zone
void Multigrid::RetrieveResult(AthenaArray<Real> &dst, int ns, int ngh)
{
  const AthenaArray<Real> &src=src_[nlev_-1];
  int is, ie, js, je, ks, ke;
  is=js=ks=0;
  ie=nx_+2*ngh_-1, je=ny_+2*ngh_-1, ke=nz_+2*ngh_-1;
#pragma ivdep
  for(int n=0; n<nvar; n++) {
    int ndst=ns+n;
#pragma ivdep
    for(int k=ngh-ngh_, mk=ks; mk<=ke; k++, mk++) {
#pragma ivdep
      for(int j=ngh-ngh_, mj=js; mj<=je; j++, mj++) {
#pragma ivdep
        for(int i=ngh-ngh_, mi=is; mi<=ie; i++, mi++)
          dst(ndst,k,j,i)=src(n,mk,mj,mi);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::ZeroClearData(int lev)
//  \brief Clear the data array with zero
void Multigrid::ZeroClearData(int lev)
{
  AthenaArray<Real> &u=u_[lev];
  std::memset(u.data(), 0, u.GetSizeInBytes());
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::ApplyPhysicalBoundaries(int lev)
//  \brief A

void Multigrid::ApplyPhysicalBoundaries(int lev)
{
  AthenaArray<Real> &dst=u_[lev];
  int ll=nlev_-1-lev;
  int ncx=nx_>>ll, ncy=ny_>>ll, ncz=nz_>>ll;
  int is=ngh_, ie=ncx_+ngh_-1, js=ngh_, je=ncy_+ngh_-1, ks=ngh_, ke=ncz_+ngh_-1;
  int bis=is-ngh_, bie=ie+ngh_, bjs=js, bje=je, bks=ks, bke=ke;
  Real dx, dy, dz;
  Real dx=rdx_/(Real)(1<<ll), dy=rdy_/(Real)(1<<ll), dz=rdz_/(Real)(1<<ll);
  Real x0=size.x1min-((Real)ngh_+0.5)*dx;
  Real y0=size.x2min-((Real)ngh_+0.5)*dy;
  Real z0=size.x3min-((Real)ngh_+0.5)*dz;
  if(MGBoundaryFunction_[INNER_X2]==NULL) bjs=js-ngh_;
  if(MGBoundaryFunction_[OUTER_X2]==NULL) bje=je+ngh_;
  if(MGBoundaryFunction_[INNER_X3]==NULL) bks=ks-ngh_;
  if(MGBoundaryFunction_[OUTER_X3]==NULL) bke=ke*ngh_;

  // Apply boundary function on inner-x1
  if (MGBoundaryFunction_[INNER_X1] != NULL)
    MGBoundaryFunction_[INNER_X1](dst, time, dt, nvar, is, ie, bjs, bje, bks, bke,
                                  x0, y0, z0, dx, dy, dz);
  // Apply boundary function on outer-x1
  if (MGBoundaryFunction_[OUTER_X1] != NULL)
    MGBoundaryFunction_[OUTER_X1](dst, time, dt, nvar, is, ie, bjs, bje, bks, bke,
                                  x0, y0, z0, dx, dy, dz);

  // Apply boundary function on inner-x2
  if (MGBoundaryFunction_[INNER_X2] != NULL)
    MGBoundaryFunction_[INNER_X2](dst, time, dt, nvar, bis, bie, js, je, bks, bke,
                                  x0, y0, z0, dx, dy, dz);
  // Apply boundary function on outer-x2
  if (MGBoundaryFunction_[OUTER_X2] != NULL)
    MGBoundaryFunction_[OUTER_X2](dst, time, dt, nvar, bis, bie, js, je, bks, bke,
                                  x0, y0, z0, dx, dy, dz);

  bjs=js-ngh_, bje=je+ngh_;
  // Apply boundary function on inner-x3
  if (MGBoundaryFunction_[INNER_X3] != NULL)
    MGBoundaryFunction_[INNER_X3](dst, time, dt, nvar, bis, bie, bjs, bje, ks, ke,
                                  x0, y0, z0, dx, dy, dz);
  // Apply boundary function on outer-x3
  if (MGBoundaryFunction_[OUTER_X3] != NULL)
    MGBoundaryFunction_[OUTER_X3](dst, time, dt, nvar, bis, bie, bjs, bje, ks, ke,
                                  x0, y0, z0, dx, dy, dz);

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::Restrict(int clev)
//  \brief Restrict the potential/density to level=clev
void Multigrid::Restrict(int clev)
{
  AthenaArray<Real> &dst=src_[clev];
  AthenaArray<Real> &src=def_[clev+1];
  int ll=nlev_-1-clev;
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+(nx_>>ll)-1, je=js+(ny_>>ll)-1, ke=ks+(nz_>>ll)-1;
#pragma ivdep
  for(int n=0; n<nvar; n++) {
#pragma ivdep
    for(int k=ks, fk=ks; k<=ke; k++, fk+=2) {
#pragma ivdep
      for(int j=js, fj=js; j<=je; j++, fj+=2) {
#pragma ivdep
        for(int i=is, fi=is; i<=ie; i++, fi+=2)
          dst(n, k, j, i)=0.125*(src(n, fk,   fj,   fi)+src(n, fk,   fj,   fi+1)
                                +src(n, fk,   fj+1, fi)+src(n, fk,   fj+1, fi+1)
                                +src(n, fk+1, fj,   fi)+src(n, fk+1, fj,   fi+1)
                                +src(n, fk+1, fj+1, fi)+src(n, fk+1, fj+1, fi+1));
      }
    }
  }
  return;
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::ProlongateAndCorrect(int clev)
//  \brief Prolongate the potential from level=clev using tri-linear interpolation
void Multigrid::ProlongateAndCorrect(int clev)
{
  const AthenaArray<Real> &src=u_[clev];
  AthenaArray<Real> &dst=def_[clev+1];
  int ll=nlev_-1-clev;
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+(nx_>>ll)-1, je=js+(ny_>>ll)-1, ke=ks+(nz_>>ll)-1;
#pragma ivdep
  for(int n=0; n<nvar; n++) {
#pragma ivdep
    for(int k=ks, fk=ks; k<=ke; k++, fk+=2) {
#pragma ivdep
      for(int j=js, fj=js; j<=je; j++, fj+=2) {
#pragma ivdep
        for(int i=is, fi=is; i<=ie; i++, fi+=2) {
          dst(n,fk  ,fj  ,fi  )+=0.015625*(27.0*src(n,k,j,i) + src(n,k-1,j-1,i-1)
                          +9.0*(src(n,k,j,i-1)+src(n,k,j-1,i)+src(n,k-1,j,i))
                          +3.0*(src(n,k-1,j-1,i)+src(n,k-1,j,i-1)+src(n,k,j-1,i-1)));
          dst(n,fk  ,fj  ,fi+1)+=0.015625*(27.0*src(n,k,j,i) + src(n,k-1,j-1,i+1)
                          +9.0*(src(n,k,j,i+1)+src(n,k,j-1,i)+src(n,k-1,j,i))
                          +3.0*(src(n,k-1,j-1,i)+src(n,k-1,j,i+1)+src(n,k,j-1,i+1)));
          dst(n,fk  ,fj+1,fi  )+=0.015625*(27.0*src(n,k,j,i) + src(n,k-1,j+1,i-1)
                          +9.0*(src(n,k,j,i-1)+src(n,k,j+1,i)+src(n,k-1,j,i))
                          +3.0*(src(n,k-1,j+1,i)+src(n,k-1,j,i-1)+src(n,k,j+1,i-1)));
          dst(n,fk+1,fj  ,fi  )+=0.015625*(27.0*src(n,k,j,i) + src(n,k+1,j-1,i-1)
                          +9.0*(src(n,k,j,i-1)+src(n,k,j-1,i)+src(n,k+1,j,i))
                          +3.0*(src(n,k+1,j-1,i)+src(n,k+1,j,i-1)+src(n,k,j-1,i-1)));
          dst(n,fk+1,fj+1,fi  )+=0.015625*(27.0*src(n,k,j,i) + src(n,k+1,j+1,i-1)
                          +9.0*(src(n,k,j,i-1)+src(n,k,j+1,i)+src(n,k+1,j,i))
                          +3.0*(src(n,k+1,j+1,i)+src(n,k+1,j,i-1)+src(n,k,j+1,i-1)));
          dst(n,fk+1,fj  ,fi+1)+=0.015625*(27.0*src(n,k,j,i) + src(n,k+1,j-1,i+1)
                          +9.0*(src(n,k,j,i+1)+src(n,k,j-1,i)+src(n,k+1,j,i))
                          +3.0*(src(n,k+1,j-1,i)+src(n,k+1,j,i+1)+src(n,k,j-1,i+1)));
          dst(n,fk  ,fj+1,fi+1)+=0.015625*(27.0*src(n,k,j,i) + src(n,k-1,j+1,i+1)
                          +9.0*(src(n,k,j,i+1)+src(n,k,j+1,i)+src(n,k-1,j,i))
                          +3.0*(src(n,k-1,j+1,i)+src(n,k-1,j,i+1)+src(n,k,j+1,i+1)));
          dst(n,fk+1,fj+1,fi+1)+=0.015625*(27.0*src(n,k,j,i) + src(n,k+1,j+1,i+1)
                          +9.0*(src(n,k,j,i+1)+src(n,k,j+1,i)+src(n,k+1,j,i))
                          +3.0*(src(n,k+1,j+1,i)+src(n,k+1,j,i+1)+src(n,k,j+1,i+1)));
        }
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::FMGProlongate(int clev)
//  \brief Prolongate the potential from level=clev for Full Multigrid cycle
void Multigrid::FMGProlongate(int clev)
{
  AthenaArray<Real> &src=u_[clev];
  AthenaArray<Real> &dst=u_[clev+1];
  int ll=nlev_-1-clev;
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+(nx_>>ll)-1, je=js+(ny_>>ll)-1, ke=ks+(nz_>>ll)-1;
#pragma ivdep
  for(int n=0; n<nvar; n++) {
#pragma ivdep
    for(int k=ks, fk=ks; k<=ke; k++, fk+=2) {
#pragma ivdep
      for(int j=js, fj=js; j<=je; j++, fj+=2) {
#pragma ivdep
        for(int i=is, fi=is; i<=ie; i++, fi+=2) {
          dst(n,fk  ,fj  ,fi  )=0.015625*(27.0*src(n,k,j,i) + src(n,k-1,j-1,i-1)
                          +9.0*(src(n,k,j,i-1)+src(n,k,j-1,i)+src(n,k-1,j,i))
                          +3.0*(src(n,k-1,j-1,i)+src(n,k-1,j,i-1)+src(n,k,j-1,i-1)));
          dst(n,fk  ,fj  ,fi+1)=0.015625*(27.0*src(n,k,j,i) + src(n,k-1,j-1,i+1)
                          +9.0*(src(n,k,j,i+1)+src(n,k,j-1,i)+src(n,k-1,j,i))
                          +3.0*(src(n,k-1,j-1,i)+src(n,k-1,j,i+1)+src(n,k,j-1,i+1)));
          dst(n,fk  ,fj+1,fi  )=0.015625*(27.0*src(n,k,j,i) + src(n,k-1,j+1,i-1)
                          +9.0*(src(n,k,j,i-1)+src(n,k,j+1,i)+src(n,k-1,j,i))
                          +3.0*(src(n,k-1,j+1,i)+src(n,k-1,j,i-1)+src(n,k,j+1,i-1)));
          dst(n,fk+1,fj  ,fi  )=0.015625*(27.0*src(n,k,j,i) + src(n,k+1,j-1,i-1)
                          +9.0*(src(n,k,j,i-1)+src(n,k,j-1,i)+src(n,k+1,j,i))
                          +3.0*(src(n,k+1,j-1,i)+src(n,k+1,j,i-1)+src(n,k,j-1,i-1)));
          dst(n,fk+1,fj+1,fi  )=0.015625*(27.0*src(n,k,j,i) + src(n,k+1,j+1,i-1)
                          +9.0*(src(n,k,j,i-1)+src(n,k,j+1,i)+src(n,k+1,j,i))
                          +3.0*(src(n,k+1,j+1,i)+src(n,k+1,j,i-1)+src(n,k,j+1,i-1)));
          dst(n,fk+1,fj  ,fi+1)=0.015625*(27.0*src(n,k,j,i) + src(n,k+1,j-1,i+1)
                          +9.0*(src(n,k,j,i+1)+src(n,k,j-1,i)+src(n,k+1,j,i))
                          +3.0*(src(n,k+1,j-1,i)+src(n,k+1,j,i+1)+src(n,k,j-1,i+1)));
          dst(n,fk  ,fj+1,fi+1)=0.015625*(27.0*src(n,k,j,i) + src(n,k-1,j+1,i+1)
                          +9.0*(src(n,k,j,i+1)+src(n,k,j+1,i)+src(n,k-1,j,i))
                          +3.0*(src(n,k-1,j+1,i)+src(n,k-1,j,i+1)+src(n,k,j+1,i+1)));
          dst(n,fk+1,fj+1,fi+1)=0.015625*(27.0*src(n,k,j,i) + src(n,k+1,j+1,i+1)
                          +9.0*(src(n,k,j,i+1)+src(n,k,j+1,i)+src(n,k+1,j,i))
                          +3.0*(src(n,k+1,j+1,i)+src(n,k+1,j,i+1)+src(n,k,j+1,i+1)));
        }
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn Real Multigrid::CalculateDefectNorm(int n, int nrm)
//  \brief calculate the residual norm

Real Multigrid::CalculateDefectNorm(int n, int nrm)
{
  AthenaArray<Real> &def=def_[nlev_-1];
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+nx_-1, je=js+ny_-1, ke=ks+nz_-1;

  Real norm=0.0;
  if(nrm==0) { // special case: max norm
#pragma ivdep
    for(int k=ks; k<=ke; k++) {
#pragma ivdep
      for(int j=js; j<=je; j++) {
#pragma ivdep
        for(int i=is; i<=ie; i++)
          norm=std::max(norm,std::fabs(def(n,k,j,i)));
      }
    }
  }
  else if (nrm==1) {
#pragma ivdep
    for(int k=ks; k<=ke; k++) {
#pragma ivdep
      for(int j=js; j<=je; j++) {
#pragma ivdep
        for(int i=is; i<=ie; i++)
          norm+=std::fabs(def(n,k,j,i));
      }
    }
  }
  else { // nrm>1 -> nrm=2
#pragma ivdep
    for(int k=ks; k<=ke; k++) {
#pragma ivdep
      for(int j=js; j<=je; j++) {
#pragma ivdep
        for(int i=is; i<=ie; i++)
          norm+=SQR(def(n,k,j,i));
      }
    }
    norm=std::sqrt(norm);
  }
  return norm;
}


//----------------------------------------------------------------------------------------
//! \fn Real Multigrid::CalculateTotalSource(int n)
//  \brief calculate the sum of the source function

Real CalculateTotalSource(int n)
{
  AthenaArray<Real> &src=src_[nlev_-1];
  Real s=0.0;
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+nx_-1, je=js+ny_-1, ke=ks+nz_-1;

#pragma ivdep
  for(int k=ks; k<=ke; k++) {
#pragma ivdep
    for(int j=js; j<=je; j++) {
#pragma ivdep
      for(int i=is; i<=ie; i++)
        s+=src(n,k,j,i);
    }
  }
  return s*rdx_*rdy_*rdz_;
}


//----------------------------------------------------------------------------------------
//! \fn Real Multigrid::SubtractMeanSource(int n, Real ave)
//  \brief subtract the mean value of the source function for periodic boundary cases

void Multigrid::SubtractMeanSource(int n, Real ave)
{
  AthenaArray<Real> &src=src_[nlev_-1];
  Real s=0.0;
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+nx_-1, je=js+ny_-1, ke=ks+nz_-1;

#pragma ivdep
  for(int k=ks; k<=ke; k++) {
#pragma ivdep
    for(int j=js; j<=je; j++) {
#pragma ivdep
      for(int i=is; i<=ie; i++)
        src(n,k,j,i)-=ave;
    }
  }
}


//----------------------------------------------------------------------------------------
//! \fn int Multigrid::GetCurrentNumberOfCells(void)
//  \brief returns the number of the cells in one direction

int Multigrid::GetCurrentNumberOfCells(void)
{
  return 1<<current_level_;
}


//----------------------------------------------------------------------------------------
//! \fn AthenaArray<Reall>& Multigrid::GetCurrentData(void)
//  \brief returns the reference of the current level

AthenaArray<Real>& Multigrid::GetCurrentData(void)
{
  return u_[current_level_];
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::Multigrid(MeshBlock *pmb, int invar, int nx, int ny, int nz,
//                                RegionSize isize, MGBoundaryFunc_t *MGBoundary)
//  \brief Multigrid constructor

void Multigrid::Multigrid(MeshBlock *pmb, int invar, int nx, int ny, int nz,
                          RegionSize isize, MGBoundaryFunc_t *MGBoundary)
{
  pmy_block_=pmb;
  size_=isize;
  nvar=invar;
  ngh_=1;
  nx_=nx, ny_=ny, nz_=nz;
  rdx_=(size_.x1max-size_.x1min)/(Real)nx;
  rdy_=(size_.x2max-size_.x2min)/(Real)ny;
  rdz_=(size_.x3max-size_.x3min)/(Real)nz;

  nlev_=0;
  int n = std::min(nx,std::min(ny, nz));
  for(int l=0; l<20; l++) {
    if((1<<l) == n) {
      nlev_=l+1;
      break;
    }
  }

  for(int i=0; i<6; i++)
    MGBoundaryFunction_[i]=MGBoundary[i];
  if(pmb!=NULL) { // not root grid
      if(pmb->pmy_mesh->multilevel==false)
        ngh_=2;
      for(int i=0; i<6; i++) {
        if(pmb->block_bcs[i]==PERIODIC_BNDRY || pmb->block_bcs[i]==BLOCK_BNDRY)
          MGBoundaryFunction_[i]=NULL;
      }
    }
  }

  // allocate arrays
  u_ = new AthenaArray<Real>[nlev_];
  src_ = new AthenaArray<Real>[nlev_];
  def_ = new AthenaArray<Real>[nlev_];
  for(int l=nlev_-1; l>=0; l++) {
    int ll=nlev_-1-l;
    int ncx=(nx>>ll)+2*ngh, ncy=(ny>>ll)+2*ngh, ncz=(nz>>ll)+2*ngh;
    u_[l].NewAthenaArray(nvar,ncz,ncy,ncx);
    src_[l].NewAthenaArray(nvar,ncz,ncy,ncx);
    def_[l].NewAthenaArray(nvar,ncz,ncy,ncx);
  }
}

//----------------------------------------------------------------------------------------
//! \fn virtual void Multigrid::~Multigrid
//  \brief Multigrid destroctor

virtual void Multigrid::~Multigrid()
{
  for(int l=0; l<nlev_; l++) {
    u_[l].DeleteAthenaArray();
    src_[l].DeleteAthenaArray();
    def_[l].DeleteAthenaArray();
  }
  delete [] u_;
  delete [] src_;
  delete [] def_;
}


//----------------------------------------------------------------------------------------
//! \fn MGPeriodicInnerX1(AthenaArray<Real> &dst,Real time, Real dt,
//                     int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
//                     Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief Periodic (default) boundary condition in the inner-X1 direction

void MGPeriodicInnerX1(AthenaArray<Real> &dst,Real time, Real dt,
                       int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
{
#pragma ivdep
  for(int n=0; n<nvar; n++) {
#pragma ivdep
    for(int k=ks, k<=ke; k++) {
#pragma ivdep
      for(int j=js; j<=je; j++) {
#pragma ivdep
        for(int i=0; i<ngh; i++)
          dst(n,k,j,is-i-1)=dst(n,k,j,ie-i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn MGPeriodicOuterX1(AthenaArray<Real> &dst,Real time, Real dt,
//                     int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
//                     Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief Periodic (default) boundary condition in the outer-X1 direction

void MGPeriodicOuterX1(AthenaArray<Real> &dst,Real time, Real dt,
                       int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
{
#pragma ivdep
  for(int n=0; n<nvar; n++) {
#pragma ivdep
    for(int k=ks, k<=ke; k++) {
#pragma ivdep
      for(int j=js; j<=je; j++) {
#pragma ivdep
        for(int i=0; i<ngh; i++)
          dst(n,k,j,ie+i+1)=dst(n,k,j,is+i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn MGPeriodicInnerX2(AthenaArray<Real> &dst,Real time, Real dt,
//                     int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
//                     Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief Periodic (default) boundary condition in the inner-X2 direction

void MGPeriodicInnerX2(AthenaArray<Real> &dst,Real time, Real dt,
                       int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
{
#pragma ivdep
  for(int n=0; n<nvar; n++) {
#pragma ivdep
    for(int k=ks, k<=ke; k++) {
#pragma ivdep
      for(int j=0; j<ngh; j++) {
#pragma ivdep
        for(int i=is; i<=ie; i++)
          dst(n,k,js-j-1,i)=dst(n,k,je-j,i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn MGPeriodicOuterX2(AthenaArray<Real> &dst,Real time, Real dt,
//                     int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
//                     Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief Periodic (default) boundary condition in the outer-X2 direction

void MGPeriodicOuterX2(AthenaArray<Real> &dst,Real time, Real dt,
                       int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
{
#pragma ivdep
  for(int n=0; n<nvar; n++) {
#pragma ivdep
    for(int k=ks, k<=ke; k++) {
#pragma ivdep
      for(int j=0; j<ngh; j++) {
#pragma ivdep
        for(int i=is; i<=ie; i++)
          dst(n,k,je+j+1,i)=dst(n,k,js+j,i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn MGPeriodicInnerX3(AthenaArray<Real> &dst,Real time, Real dt,
//                     int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
//                     Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief Periodic (default) boundary condition in the inner-X3 direction

void MGPeriodicInnerX3(AthenaArray<Real> &dst, Real time, Real dt,
                       int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
{
#pragma ivdep
  for(int n=0; n<nvar; n++) {
#pragma ivdep
    for(int k=0, k<=ngh; k++) {
#pragma ivdep
      for(int j=js; j<=je; j++) {
#pragma ivdep
        for(int i=is; i<=ie; i++)
          dst(n,ks-k-1,j,i)=dst(n,ke-k,j,i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn MGPeriodicOuterX3(AthenaArray<Real> &dst,Real time, Real dt,
//                int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
//                Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief Periodic (default) boundary condition in the outer-X3 direction

void MGPeriodicOuterX3(AthenaArray<Real> &dst,Real time, Real dt,
                       int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
{
#pragma ivdep
  for(int n=0; n<nvar; n++) {
#pragma ivdep
    for(int k=0, k<=ngh; k++) {
#pragma ivdep
      for(int j=js; j<=je; j++) {
#pragma ivdep
        for(int i=is; i<=ie; i++)
          dst(n,ke+k+1,j,i)=dst(n,ks+k,j,i);
      }
    }
  }
  return;
}

