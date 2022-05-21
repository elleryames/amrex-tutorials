
#include <AMReX.H>
#include "MyTest.H"

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    {
        BL_PROFILE("main");
        MyTest mytest;
        mytest.newtonSolve();
        mytest.writePlotfile();
    }

    amrex::Finalize();
}
