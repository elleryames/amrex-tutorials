#include "MyTest.H"

#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;

MyTest::MyTest ()
{
    readParameters();
    initData();
}


void
MyTest::newtonSolve ()
{
    // Initial set up for solver
    LPInfo info;
    info.setAgglomeration(agglomeration);
    info.setConsolidation(consolidation);
    info.setSemicoarsening(semicoarsening);
    info.setMaxCoarseningLevel(max_coarsening_level);
    info.setMaxSemicoarseningLevel(max_semicoarsening_level);

    // linear solver params
    const Real tol_rel = 1.e-10;
    const Real tol_abs = 0.0;

    // Newton-Raphson iteration
    Real nr_res = 1.0;

    const int nlevels = geom.size();

    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        amrex::Print() << "Level " << ilev << std::endl;

        MLABecLaplacian mlabec({geom[ilev]}, {grids[ilev]}, {dmap[ilev]}, info);

        mlabec.setMaxOrder(linop_maxorder);

        // This is a 2d problem with homogeneous Dirichlet BC
        mlabec.setDomainBC({LinOpBCType::Dirichlet,
                                LinOpBCType::Dirichlet},
                            {LinOpBCType::Dirichlet,
                                LinOpBCType::Dirichlet});

        if (ilev > 0) {
            mlabec.setCoarseFineBC(&solution[ilev-1], ref_ratio);
        }

        // for problem with pure homogeneous Dirichlet BC, we could pass a nullptr
        mlabec.setLevelBC(0, nullptr);

        mlabec.setScalars(ascalar, bscalar);

        Array<MultiFab,AMREX_SPACEDIM> face_bcoef;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            const BoxArray& ba = amrex::convert(bcoef[ilev].boxArray(),
                                                IntVect::TheDimensionVector(idim));
            face_bcoef[idim].define(ba, bcoef[ilev].DistributionMap(), 1, 0);
        }
        amrex::average_cellcenter_to_face(GetArrOfPtrs(face_bcoef),
                                            bcoef[ilev], geom[ilev]);
        mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(face_bcoef));

        // write in terms of for with maxiter
        for (int nriter = 0; nriter < nr_maxiter; nriter++)
        {
            // Update linear problem
            // Note that for this example, the alpha coefficient is the solution from the previous iterate.
            mlabec.setACoeffs(0, solution[ilev]);
    
            // create multigrid solver for updated equation
            MLMG mlmg(mlabec);
            mlmg.setMaxIter(max_iter);
            mlmg.setMaxFmgIter(max_fmg_iter);
            mlmg.setVerbose(verbose);
            mlmg.setBottomVerbose(bottom_verbose);

            // compute residual
            mlmg.compResidual ({&residual[ilev]}, {&solution[ilev]}, 
                                {&rhs[ilev]});
            auto dx = geom[ilev].CellSize();
            nr_res = residual[ilev].norm2(0) * dx[0] * dx[1];
            amrex::Print() << " Newton step: " << nriter
                           << " residual: " << nr_res << std::endl;

            // Compute norm of residual and assign to nr_res
            // nr_res = *(std::max_element(nr_err.begin(), nr_err.end()));
            if (nr_res < nr_tol) {break;}

            // solve linear problem
            mlmg.solve({&linsol[ilev]}, {&rhs[ilev]}, tol_rel, tol_abs);
                
            // Update solution
            MultiFab::Add(solution[ilev], linsol[ilev], 0, 0, 1, 0);
        }
    }
}

// void
// MyTest::linearSolve (MLLinOp& mlabec)
// {
//     // Do we need to create a new MLMG solver each iteration or can we reuse one?
//     // MLMG mlmg(mlabec);
//     // mlmg.setMaxIter(max_iter);
//     // mlmg.setMaxFmgIter(max_fmg_iter);
//     // mlmg.setVerbose(verbose);
//     // mlmg.setBottomVerbose(bottom_verbose);

//     // mlmg.solve({&solution[ilev]}, {&rhs[ilev]}, tol_rel, tol_abs);
// }

// Real 
// MyTest::computeResidual ()
// {
//     // Do I have access to the laplacian?
//     // Can objects be subtracted directly, or do I have to compute this pointwise? 
//     // How to call an L2 norm? 
// }

void
MyTest::readParameters ()
{
    ParmParse pp;
    pp.query("max_level", max_level);
    pp.query("ref_ratio", ref_ratio);
    pp.query("n_cell", n_cell);
    pp.query("max_grid_size", max_grid_size);

    pp.query("composite_solve", composite_solve);

    pp.query("verbose", verbose);
    pp.query("bottom_verbose", bottom_verbose);
    pp.query("max_iter", max_iter);
    pp.query("max_fmg_iter", max_fmg_iter);
    pp.query("linop_maxorder", linop_maxorder);
    pp.query("agglomeration", agglomeration);
    pp.query("consolidation", consolidation);
    pp.query("semicoarsening", semicoarsening);
    pp.query("max_coarsening_level", max_coarsening_level);
    pp.query("max_semicoarsening_level", max_semicoarsening_level);

    pp.query("nr_tol", nr_tol);
    pp.query("nr_maxiter", nr_maxiter);
}

void
MyTest::initData ()
{
    int nlevels = max_level + 1;
    geom.resize(nlevels);
    grids.resize(nlevels);
    dmap.resize(nlevels);

    solution.resize(nlevels);
    rhs.resize(nlevels);
    linsol.resize(nlevels);
    residual.resize(nlevels);
    acoef.resize(nlevels);
    bcoef.resize(nlevels);

    // Set up Computational Domain
    // Domain is box [0,1]x[0,1]
    RealBox rb({0.,0.}, {1.,1.});
    // Non-periodic BC
    Array<int,AMREX_SPACEDIM> is_periodic{0,0};
    // Cartesian coords
    int coordsys = 0;
    Geometry::Setup(&rb, coordsys, is_periodic.data());
    Box domain0(IntVect{0,0}, IntVect{n_cell-1,n_cell-1});
    Box domain = domain0;
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        geom[ilev].define(domain);
        domain.refine(ref_ratio);
    }

    domain = domain0;
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        grids[ilev].define(domain);
        grids[ilev].maxSize(max_grid_size);
        domain.grow(-n_cell/4);   // fine level cover the middle of the coarse domain
        domain.refine(ref_ratio);
    }

    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        dmap[ilev].define(grids[ilev]);
        // mf.define (BoxArray&, DistributionMapping&, int numvar, IntVect& ngrow)
        solution      [ilev].define(grids[ilev], dmap[ilev], 1, 1);
        rhs           [ilev].define(grids[ilev], dmap[ilev], 1, 0);
        linsol        [ilev].define(grids[ilev], dmap[ilev], 1, 1);
        residual      [ilev].define(grids[ilev], dmap[ilev], 1, 0);
        if (!acoef.empty()) {
            acoef[ilev].define(grids[ilev], dmap[ilev], 1, 0);
            bcoef[ilev].define(grids[ilev], dmap[ilev], 1, 1);
        }
    }

    // Initialize problem
    initProbABecLaplacian();
    
}

