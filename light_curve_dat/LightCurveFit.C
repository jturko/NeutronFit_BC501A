
TF1 * LightCurveFit(int i)
{

    TGraph * graph = new TGraph(Form("v%d.dat",i));
    TF1 * fitfunc = new TF1("fitfunc","[0]*x - [1]*(1.-TMath::Exp(-[2]*TMath::Power(x,[3])))",0.001,10);
    //fitfunc->SetParameters( 0.7869, 3.1755, 0.2069, 0.9558 );
    fitfunc->SetParameters(0.9393, 4.7507, 0.1722 , 1.0830);
    fitfunc->SetParNames("a1","a2","a3","a4");
    graph->Fit(fitfunc);

    //fitfunc->Draw();
    graph->Draw("A*");
    
    return fitfunc;

}
