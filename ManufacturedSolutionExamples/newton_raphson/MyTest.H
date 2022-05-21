#ifndef MY_TEST_H_
#define MY_TEST_H_

#include <AMReX_MLMG.H>

#ifdef AMREX_USE_HYPRE
#include <AMReX_Hypre.H>
#endif

class MyTest
{
public:

    MyTest ();

    void newtonSolve ();
    void writePlotfile () const;

public: // make these public for cuda
    void initProbABecLaplacian ();

private:

    void readParameters ();
    void initData ();
    // void linearSolve (MLLinOp& mlabec);
    // amrex::Real computeResidual ();

    int max_level = 1;
    int ref_ratio = 2;
    int n_cell = 128;
    int max_grid_size = 64;

    bool composite_solve = false;

    // For MLMG solver
    int verbose = 2;
    int bottom_verbose = 0;
    int max_iter = 100;
    int max_fmg_iter = 0;
    int linop_maxorder = 2;
    bool agglomeration = true;
    bool consolidation = true;
    bool semicoarsening = false;
    int max_coarsening_level = 30;
    int max_semicoarsening_level = 0;
    bool use_hypre = false;
    bool use_petsc = false;

    // For Newton solver
    amrex::Real nr_tol = 1e-6;
    int  nr_maxiter = 2;

    amrex::Vector<amrex::Geometry> geom;
    amrex::Vector<amrex::BoxArray> grids;
    amrex::Vector<amrex::DistributionMapping> dmap;

    amrex::Vector<amrex::MultiFab> solution;
    amrex::Vector<amrex::MultiFab> linsol; // for solution to linear problem
    amrex::Vector<amrex::MultiFab> residual; // residual of nonlinear problem
    amrex::Vector<amrex::MultiFab> rhs;
    amrex::Vector<amrex::MultiFab> acoef;
    amrex::Vector<amrex::MultiFab> bcoef;

    amrex::Real ascalar = -2.0;
    amrex::Real bscalar = -1.0;
};

#endif
