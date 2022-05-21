#ifndef INIT_PROB_K_H_
#define INIT_PROB_K_H_

#include <AMReX_FArrayBox.H>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void actual_init_bcoef (int i, int j, int k,
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

    constexpr amrex::Real s2 = 0.25;
    amrex::Real recs4 = 1./s2/s2;

    amrex::Real r = prob_lo[0] + dx[0] * (i + 0.5);
    amrex::Real z = prob_lo[1] + dx[1] * (j + 0.5);

    // Robin BC applied on outer (non-axis) faces. 
    bool robin_cell = (i > dhi.x) || (j > dhi.y) || (j < dlo.y);
    if (robin_cell) {
        robin_a(i,j,k) = 1.0;
        robin_b(i,j,k) = 1.0; 
        amrex::Real dphidr = -2. * r / s2 * exact(i,j,k);
        amrex::Real dphidz = -2. * z / s2 * exact(i,j,k);
        amrex::Real dphidn = r * dphidr + z * dphidz;
        robin_f(i,j,k) = robin_a(i,j,k) * exact(i,j,k) 
                            + robin_b(i,j,k) * dphidn;
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void actual_init_abeclap (int i, int j, int k,
                          amrex::Array4<amrex::Real      > const& rhs,
                          amrex::Array4<amrex::Real      > const& exact,
                          amrex::Array4<amrex::Real      > const& alpha,
                          amrex::Real a, amrex::Real b,
                          amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_lo,
                          amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_hi,
                          amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
                          const amrex::Dim3& dlo, 
                          const amrex::Dim3& dhi)
{
    constexpr amrex::Real pi = 3.1415926535897932;
    constexpr amrex::Real s2 = 0.25;
    amrex::Real recs4 = 1./s2/s2;

    amrex::Real r = prob_lo[0] + dx[0] * (i + 0.5);
    amrex::Real z = prob_lo[1] + dx[1] * (j + 0.5);

    alpha(i,j,k) = 0.0;
    exact(i,j,k) = std::exp(-(r*r + z*z)/s2);
    rhs(i,j,k) = (6. * s2 - 4. * (r*r + z*z) ) * recs4 * exact(i,j,k);

    // // Robin BC applied on outer (non-axis) faces. 
    // bool robin_cell = (i > dhi.x) || (j > dhi.y) || (j < dlo.y);
    // if (robin_cell) {
    //     robin_a(i,j,k) = 1.0;
    //     robin_b(i,j,k) = 1.0; 
    //     amrex::Real dphidr = -2. * r / s2 * exact(i,j,k);
    //     amrex::Real dphidz = -2. * z / s2 * exact(i,j,k);
    //     amrex::Real dphidn = r * dphidr + z * dphidz;
    //     robin_f(i,j,k) = robin_a(i,j,k) * exact(i,j,k) 
    //                         + robin_b(i,j,k) * dphidn;
    // }
}

#endif