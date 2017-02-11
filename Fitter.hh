
#ifndef FITTER_H
#define FITTER_H
#endif

#include "NeutronFit_BC501A.cc"
#include "vec.hh"

#include "Math/GSLMinimizer.h"
#include "Math/Functor.h"

class Fitter 
{

public:
    Fitter();
    ~Fitter();

    Fitter(int a);
    Fitter(int a, int b);
    Fitter(int a, int b, int c);
    Fitter(int a, int b, int c, int d);
    Fitter(int a, int b, int c, int d, int e, int f);
    Fitter(int a, int b, int c, int d, int e, int f, int g, int h);
    Fitter(int a, int b, int c, int d, int e, int f, int g, int h, int i, int j);

    void Draw();
    void Run(double a1=0.639, double a2=1.462, double a3=0.373, double a4=0.968, double carbon=0);

    bool Check(int i) { if(i<=-1||i>=GetNumberOfNeutronFit_BC501As()) return false; else return true; }
    
    NeutronFit_BC501A * GetNeutronFit_BC501A(int i) { if(Check(i)) return &fNeutronFit_BC501AVector[i]; else{ std::cout << "out of bounds!" << std::endl; return NULL;} }
    
    void SetNextNeutronFit_BC501A(NeutronFit_BC501A hfit) { fNeutronFit_BC501AVector.push_back(hfit); fRunNumVector.push_back(hfit.GetRunNum());  }
    void SetNextNeutronFit_BC501A(int i) {
        NeutronFit_BC501A * hfit = new NeutronFit_BC501A(i);
        SetNextNeutronFit_BC501A(*hfit);
    }    

    int GetNumberOfNeutronFit_BC501As() { return fNeutronFit_BC501AVector.size(); }
    
    void DrawNeutronFit_BC501A(int i) { 
        if(!Check(i)) { std::cout << "out of bounds!" << std::endl; return; }
        if(fCanvas) delete fCanvas; 
        fNeutronFit_BC501AVector.at(i).Draw(); 
    }
    
    void SortRun(int num) { 
        if(Check(num)) { 
            fNeutronFit_BC501AVector.at(num).Sort(fParameters); 
        }
        else{ std::cout << "out of bounds!" << std::endl;} 
    }
    
    void SortAllRuns() { 
        PrintParameters();
        for(int num=0; num<GetNumberOfNeutronFit_BC501As(); num++) SortRun(num); 
        DoChi2();
    }
    
    void Print() { for(int num=0; num<GetNumberOfNeutronFit_BC501As(); num++) std::cout << "Run# = " << fRunNumVector.at(num) << " ; Energy = " << fNeutronFit_BC501AVector.at(num).GetEnergy() << std::endl; } 

    void SetParameters(double a1, double a2, double a3, double a4, double carbon) {
        fParameters[0]=a1;
        fParameters[1]=a2;
        fParameters[2]=a3;
        fParameters[3]=a4;
        fParameters[4]=carbon;
        for(int i=0; i<GetNumberOfNeutronFit_BC501As(); i++) fNeutronFit_BC501AVector.at(i).SetParameters(fParameters);
    }
    void SetParameters(double * par) { // expects a par array w/ 5 elements
        for(int i=0; i<5; i++) fParameters[i] = par[i];
        for(int i=0; i<GetNumberOfNeutronFit_BC501As(); i++) fNeutronFit_BC501AVector.at(i).SetParameters(fParameters);
    }
    void SetOffset(double offset) {
        for(int i=0; i<GetNumberOfNeutronFit_BC501As(); i++) fNeutronFit_BC501AVector.at(i).SetOffset(offset);
    }    
    void SetSmearingCoeff(double A, double B, double C) {
        for(int i=0; i<GetNumberOfNeutronFit_BC501As(); i++) fNeutronFit_BC501AVector.at(i).SetSmearingCoeff(A,B,C);
    }
    double GetOffset() { return fNeutronFit_BC501AVector.at(0).GetOffset(); }
    double GetSmearingCoeff(int i) { return fNeutronFit_BC501AVector.at(0).GetSmearingCoeff(i); }

    void PrintParameters() { 
        std::cout << "     a1 = " << fParameters[0] << std::endl;
        std::cout << "     a2 = " << fParameters[1] << std::endl;
        std::cout << "     a3 = " << fParameters[2] << std::endl;
        std::cout << "     a4 = " << fParameters[3] << std::endl;
        std::cout << " carbon = " << fParameters[4] << std::endl;
        std::cout << "      A = " << fNeutronFit_BC501AVector.at(0).GetSmearingCoeff(0) << std::endl;
        std::cout << "      B = " << fNeutronFit_BC501AVector.at(0).GetSmearingCoeff(1) << std::endl;
        std::cout << "      C = " << fNeutronFit_BC501AVector.at(0).GetSmearingCoeff(2) << std::endl;
        std::cout << " offset = " << fNeutronFit_BC501AVector.at(0).GetOffset() << std::endl;
    }

    void InitializeParameters();

    //void NelderMead(int itermax = 50);
    //vec NelderMead(vec initial_vec, int itermax = 50);
    //void NelderMead3(double a1=0.639, double a2=1.462, double a3=0.373, double a4=0.968, double carbon=0, int itermax=50);    
    vec NelderMead(double a1=0.639, double a2=1.462, double a3=0.373, double a4=0.968, double carbon=0, int itermax=50);
    vec NelderMead(vec input, int itermax=50);
    
    double DoChi2() { 
        fSum = 0.;
        fSum2 = 0.;
        for(int i=0;i<GetNumberOfNeutronFit_BC501As();i++) fSum += fNeutronFit_BC501AVector.at(i).DoChi2();
        for(int i=0;i<GetNumberOfNeutronFit_BC501As();i++) fSum2 += fNeutronFit_BC501AVector.at(i).DoChi2() * fNeutronFit_BC501AVector.at(i).DoChi2();
        fSum /= double(GetNumberOfNeutronFit_BC501As());
        fSum2 /= double(GetNumberOfNeutronFit_BC501As());    
    
        return fSum2; // CHANGE THIS TO SET WHAT WE WILL MINIMIZE [ chi2 or (chi2)^2 ]
    }
        
    double nm_val(double * par) {
        SetParameters(par);
        SortAllRuns();
        return DoChi2();
    }
    double nm_val(vec v) {
        SetParameters(v.par_array());
        SortAllRuns();
        return DoChi2();
    }
    
    double FitValue(const double * par) {
        double mypar[5];
        for(int i=0; i<5; i++) mypar[i]=par[i];
        SetParameters(mypar);
        //if(DidParametersChange(mypar)) SortAllRuns();
        SortAllRuns();
        return DoChi2();
    }
    bool DidParametersChange(double * par) {
        for(int i=0; i<5; i++) {
            if(TMath::Abs(fParameters[i] - par[i] > 0.00001)) return true;
        }
        return false;
    }    
    

    void DrawToFile(std::string name);

    int Minimize();

    std::vector<NeutronFit_BC501A> fNeutronFit_BC501AVector;   
    std::vector<int> fRunNumVector;

    double fParameters[5];   
 
    TCanvas * fCanvas;
    
    double fSum;
    double fSum2;


};

