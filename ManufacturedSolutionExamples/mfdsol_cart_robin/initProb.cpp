
#include "MyTest.H"
#include "initProb_K.H"

using namespace amrex;

void
MyTest::initProbABecLaplacian ()
{
    for (int ilev = 0; ilev <= max_level; ++ilev)
    {
        const auto prob_lo = geom[ilev].ProbLoArray();
        const auto prob_hi = geom[ilev].ProbHiArray();
        const auto dx      = geom[ilev].CellSizeArray();
        const auto dlo     = amrex::lbound(geom[ilev].Domain());
        const auto dhi     = amrex::ubound(geom[ilev].Domain());
        auto a = ascalar;
        auto b = bscalar;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(rhs[ilev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const Box& gbx = mfi.growntilebox(1);

            auto rhsfab = rhs[ilev].array(mfi);
            auto exactfab = exact_solution[ilev].array(mfi);
            auto acoeffab = acoef[ilev].array(mfi);
            auto bcoeffab = bcoef[ilev].array(mfi);
            const auto& rafab  = robinbc_a[ilev].array(mfi);
            const auto& rbfab  = robinbc_b[ilev].array(mfi);
            const auto& rffab  = robinbc_f[ilev].array(mfi);

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                actual_init_bx(i,j,k,
                                rhsfab,exactfab,acoeffab,
                                prob_lo,prob_hi,dx);
            });

            amrex::ParallelFor(gbx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                actual_init_gbx(i,j,k,
                                bcoeffab,exactfab,
                                rafab,rbfab,rffab,
                                prob_lo,prob_hi,dx,dlo,dhi);
            });
        }

        solution[ilev].setVal(0.0);
    }
}