#ifndef INIT_PROB_K_H_
#define INIT_PROB_K_H_

#include <AMReX_FArrayBox.H>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void actual_init_bcoef (int i, int j, int k,
                        amrex::Array4<amrex::Real> const& beta,
                        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_lo,
                        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_hi,
                        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx)
{
    beta(i,j,k) = 1.0;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void actual_init_abeclap (int i, int j, int k,
                          amrex::Array4<amrex::Real      > const& rhs,
                          amrex::Array4<amrex::Real      > const& alpha,
                          amrex::Real a, amrex::Real b,
                          amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_lo,
                          amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_hi,
                          amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx)
{
    amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
    amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);

    alpha(i,j,k) = 0.;
    rhs(i,j,k) = -x*(1.-x) * y*(1.-y);
}

#endif