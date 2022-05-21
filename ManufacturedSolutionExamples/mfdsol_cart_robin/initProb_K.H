#ifndef INIT_PROB_K_H_
#define INIT_PROB_K_H_

#include <AMReX_FArrayBox.H>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void actual_init_gbx (int i, int j, int k,
                        amrex::Array4<amrex::Real > const& beta,
                        amrex::Array4<amrex::Real   const > const& exact,
                        amrex::Array4<amrex::Real > const& robin_a,
                        amrex::Array4<amrex::Real > const& robin_b,
                        amrex::Array4<amrex::Real > const& robin_f,
                        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_lo,
                        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_hi,
                        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,                         
                        const amrex::Dim3& dlo, 
                        const amrex::Dim3& dhi)
{
    beta(i,j,k) = 1.0;

    constexpr amrex::Real pi = 3.1415926535897932;

    amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
    amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);

    // Robin BC applied on x=1 face. 
    bool robin_cell = (i > dhi.x);
    if (robin_cell) {
        robin_a(i,j,k) = 1.0;
        robin_b(i,j,k) = 1.0; 
        amrex::Real dphidn = std::cos(x)*std::sin(pi*y);
        robin_f(i,j,k) = robin_a(i,j,k) * exact(i,j,k) 
                            + robin_b(i,j,k) * dphidn;
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void actual_init_bx (int i, int j, int k,
                          amrex::Array4<amrex::Real      > const& rhs,
                          amrex::Array4<amrex::Real      > const& exact,
                          amrex::Array4<amrex::Real      > const& alpha,
                          amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_lo,
                          amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_hi,
                          amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx)
{
    constexpr amrex::Real pi = 3.1415926535897932;

    amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
    amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);

    alpha(i,j,k) = 0.0;
    exact(i,j,k) = std::sin(x)*std::sin(pi*y);
    rhs(i,j,k) = (1.+ pi*pi)*exact(i,j,k);

    // // Robin BC applied on x=1 face. 
    // bool robin_cell = (i > dhi.x);
    // if (robin_cell) {
    //     robin_a(i,j,k) = 1.0;
    //     robin_b(i,j,k) = 1.0; 
    //     amrex::Real dphidn = std::cos(x)*std::sin(pi*y);
    //     robin_f(i,j,k) = robin_a(i,j,k) * exact(i,j,k) 
    //                         + robin_b(i,j,k) * dphidn;
    // }
}

#endif