
#include "Fitter.cc"

int main()
{
    double energy_vector[] =
    {
        0.66201 ,0.80879    ,1.05264    ,1.21767    ,2.96445    ,2.86088    ,2.69813    ,2.48962    ,2.25152    ,
        2.00070 ,1.75268    ,1.52003    ,2.27882    ,2.60832    ,2.95748    ,3.30884    ,3.64108    ,3.93118    ,
        4.15713 ,4.30074    ,1.98127    ,1.72261    ,1.50509    ,1.50509    ,1.72261    ,1.98127    ,0.81657    ,
        0.95882 ,1.11930    ,1.29208    ,1.46825    ,1.63655    ,1.78465    ,1.90062    ,1.97458    ,0.69544    ,
        0.59575 ,0.51611    ,0.45421    ,0.38915    ,0.68765    ,0.65189    ,0.59641    ,0.48874    ,0.37184    ,
        0.26609 ,0.18371    ,0.12684    ,0.10089    ,0.07506    ,0.05996    ,0.05996    ,7.90894    ,7.64353    ,
        7.22603 ,6.69019    ,6.07682    ,5.42853    ,4.78473    ,4.17760    ,3.62976    ,3.15389    ,2.75400    ,
        2.42779
    };

    int run_vector[] =
    {
        50  ,49  ,48  ,47  ,46  ,45  ,44  ,39  ,38  ,43  ,37  ,36  ,42  ,41  ,0   ,40  ,35  ,1   ,26  ,27  ,2   ,
        28  ,3   ,29  ,30  ,22  ,23  ,11  ,31  ,21  ,24  ,10  ,32  ,33  ,34  ,20  ,25  ,9   ,8   ,12  ,63  ,7   ,
        13  ,6   ,62  ,5   ,14  ,4   ,61  ,15  ,60  ,16  ,17  ,18  ,59  ,19  ,58  ,57  ,56  ,55  ,54  ,53  ,52  
    };

    std::ofstream outfile;
    outfile.open("outfile.csv");   

    //vec v(0.09375 , 0.0751694 , 0.013922 , 0.749069 , 3.2 , 0.170967 , 0.952316 , 0.0001);
    //for(int i=62; i>=0; i--) {
    vec v;
    for(int i=62; i>=0; i--) {
        Fitter fit(run_vector[i]);
        fit.SetSmearingCoeff(0.203312 , 0.148887 , 0.000109966);
        vec v_new;
        if(i==62) v_new = fit.NelderMead();
        else      v_new = fit.NelderMead(v,50);
        v = v_new;
        outfile << run_vector[i] << " , ";
        outfile << energy_vector[run_vector[i]] << " , "; 
        for(int i=0; i<v.size(); i++) outfile << v.at(i) << " , ";
        outfile << std::endl;
    }
    
    outfile.close();

    return 0;
}
