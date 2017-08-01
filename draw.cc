
#include "vec.hh"
#include "NeutronFit_BC501A.hh"
#include "ProtonFitter.hh"

int main()
{
    ProtonFitter fit;
 
    fit.LoadAll(true);
    //for(int i=0; i<51; i++) fit.SetNextNeutronFit_BC501A(i);
    //for(int i=52; i<64; i++) fit.SetNextNeutronFit_BC501A(i);

    fit.SetSimSortMax(-1);    
    fit.SetChi2Method(2);
    
    //fit.Print();

    //fit.SetParameters(0.81921 , 3.32790 , 0.21158 , 0.96077 , 0.00635 , 0.21673 , 0.10344 , 0.00363);
    //fit.SetParameters(0.81921 , 3.32790 , 0.21158 , 0.96077 , 0 , 0.21673 , 0.10344 , 0.00363);
    //fit.SetOffset(10.45127);

    // output from fit edge script
    //fit.SetParameters(4.84081e-01 , 5.41148e-01 , 8.50818e-01 , 1.16615e+00 , 0 , 0 , 0 , 0);
    fit.SetParameters(4.84081e-01 , 5.41148e-01 , 8.50818e-01 , 1.16615e+00 , 0 , 0.38586 , 0.05694 , 0.00001);
    fit.SetOffset(0);

    fit.SortAllRunsMT();
    
    fit.PrintChi2();

    fit.DrawToFile("draw.pdf");

    return 0;

}
