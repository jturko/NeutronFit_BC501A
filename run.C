
#include "Fitter.cc"

Fitter run()
{
    Fitter fit;
    //int v[] = { 50,48,37,2,9,4,17,58,55,52 };
    //for(int i=0; i<10; i++) fit.SetNextNeutronFit_BC501A(v[i]);
    for(int i=0; i<64; i++) fit.SetNextNeutronFit_BC501A(i);
    fit.SetSmearingCoeff(0.123,0.125,0.0074);
    fit.SetParameters(0.646189 , 1.46985 , 0.375871 , 0.967337, 0);
    fit.PrintParameters();
    fit.SortAllRuns();
    fit.DrawToFile("run_out.pdf");
    return fit;
}

