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
MyTest::solve ()
{
    LPInfo info;
    info.setAgglomeration(agglomeration);
    info.setConsolidation(consolidation);
    info.setSemicoarsening(semicoarsening);
    info.setMaxCoarseningLevel(max_coarsening_level);
    info.setMaxSemicoarseningLevel(max_semicoarsening_level);

    const Real tol_rel = 1.e-10;
    const Real tol_abs = 0.0;

    const int nlevels = geom.size();

    if (composite_solve)
    {

        MLABecLaplacian mlabec(geom, grids, dmap, info);

        mlabec.setMaxOrder(linop_maxorder);

        // This is a 3d problem with Robin BC
        mlabec.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                        LinOpBCType::Dirichlet,
                                        LinOpBCType::Dirichlet)},
                            {AMREX_D_DECL(LinOpBCType::Robin,
                                        LinOpBCType::Dirichlet,
                                        LinOpBCType::Dirichlet)});                                     

        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            mlabec.setLevelBC(ilev, {&solution[ilev]}, 
                                    {&robinbc_a[ilev]}, 
                                    {&robinbc_b[ilev]}, 
                                    {&robinbc_f[ilev]});
        }

        mlabec.setScalars(ascalar, bscalar);

        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            mlabec.setACoeffs(ilev, acoef[ilev]);

            Array<MultiFab,AMREX_SPACEDIM> face_bcoef;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                const BoxArray& ba = amrex::convert(bcoef[ilev].boxArray(),
                                                    IntVect::TheDimensionVector(idim));
                face_bcoef[idim].define(ba, bcoef[ilev].DistributionMap(), 1, 0);
            }
            amrex::average_cellcenter_to_face(GetArrOfPtrs(face_bcoef),
                                              bcoef[ilev], geom[ilev]);
            mlabec.setBCoeffs(ilev, amrex::GetArrOfConstPtrs(face_bcoef));
        }

        MLMG mlmg(mlabec);
        mlmg.setMaxIter(max_iter);
        mlmg.setMaxFmgIter(max_fmg_iter);
        mlmg.setVerbose(verbose);
        mlmg.setBottomVerbose(bottom_verbose);
#ifdef AMREX_USE_HYPRE
        if (use_hypre) {
            mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
            mlmg.setHypreInterface(hypre_interface);
        }
#endif
#ifdef AMREX_USE_PETSC
        if (use_petsc) {
            mlmg.setBottomSolver(MLMG::BottomSolver::petsc);
        }
#endif

        mlmg.solve(GetVecOfPtrs(solution), GetVecOfConstPtrs(rhs), tol_rel, tol_abs);
    }
    else
    {
        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            MLABecLaplacian mlabec({geom[ilev]}, {grids[ilev]}, {dmap[ilev]}, info);

            mlabec.setMaxOrder(linop_maxorder);

            // This is a 3d problem with Robin BC
            mlabec.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                             LinOpBCType::Dirichlet,
                                             LinOpBCType::Dirichlet)},
                               {AMREX_D_DECL(LinOpBCType::Robin,
                                             LinOpBCType::Dirichlet,
                                             LinOpBCType::Dirichlet)});

            if (ilev > 0) {
                mlabec.setCoarseFineBC(&solution[ilev-1], ref_ratio);
            }

            mlabec.setLevelBC(0, {&solution[ilev]}, 
                                    {&robinbc_a[ilev]}, 
                                    {&robinbc_b[ilev]}, 
                                    {&robinbc_f[ilev]});

            mlabec.setScalars(ascalar, bscalar);

            mlabec.setACoeffs(0, acoef[ilev]);

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

            MLMG mlmg(mlabec);
            mlmg.setMaxIter(max_iter);
            mlmg.setMaxFmgIter(max_fmg_iter);
            mlmg.setVerbose(verbose);
            mlmg.setBottomVerbose(bottom_verbose);
#ifdef AMREX_USE_HYPRE
            if (use_hypre) {
                mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
                mlmg.setHypreInterface(hypre_interface);
            }
#endif
#ifdef AMREX_USE_PETSC
            if (use_petsc) {
                mlmg.setBottomSolver(MLMG::BottomSolver::petsc);
            }
#endif

            mlmg.solve({&solution[ilev]}, {&rhs[ilev]}, tol_rel, tol_abs);
        }
    }
}

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

#ifdef AMREX_USE_HYPRE
    pp.query("use_hypre", use_hypre);
    pp.query("hypre_interface", hypre_interface_i);
    if (hypre_interface_i == 1) {
        hypre_interface = Hypre::Interface::structed;
    } else if (hypre_interface_i == 2) {
        hypre_interface = Hypre::Interface::semi_structed;
    } else {
        hypre_interface = Hypre::Interface::ij;
    }
#endif
#ifdef AMREX_USE_PETSC
    pp.query("use_petsc", use_petsc);
#endif
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!(use_hypre && use_petsc),
                                     "use_hypre & use_petsc cannot be both true");
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
    exact_solution.resize(nlevels);
    acoef.resize(nlevels);
    bcoef.resize(nlevels);
    robinbc_a.resize(nlevels);
    robinbc_b.resize(nlevels);
    robinbc_f.resize(nlevels);

    // Set up Computational Domain
    // Domain is box [0,1]x[0,1]
    RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.,1.,1.)});
    // Non-periodic BC
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};
    // Cartesian coords
    int coordsys = 0;
    Geometry::Setup(&rb, coordsys, is_periodic.data());
    Box domain0(IntVect{AMREX_D_DECL(0,0,0)}, 
                IntVect{AMREX_D_DECL(n_cell-1,n_cell-1,n_cell-1)});
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
        solution      [ilev].define(grids[ilev], dmap[ilev], 1, 1);
        rhs           [ilev].define(grids[ilev], dmap[ilev], 1, 0);
        exact_solution[ilev].define(grids[ilev], dmap[ilev], 1, 0);
        robinbc_a     [ilev].define(grids[ilev], dmap[ilev], 1, 1);
        robinbc_b     [ilev].define(grids[ilev], dmap[ilev], 1, 1);
        robinbc_f     [ilev].define(grids[ilev], dmap[ilev], 1, 1);
        if (!acoef.empty()) {
            acoef[ilev].define(grids[ilev], dmap[ilev], 1, 0);
            bcoef[ilev].define(grids[ilev], dmap[ilev], 1, 1);
        }
    }

    // Initialize problem
    initProbABecLaplacian();
    
}

