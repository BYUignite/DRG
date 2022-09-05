#include "psr.h"

#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;
using namespace Cantera;

////////////////////////////////////////////////////////////////////////////////

int main() {
        
    //--------------- initialize cantera

    //auto sol = newSolution("H2O2.yaml", "", "None");
    auto sol = newSolution("gri30.yaml", "", "None");
    auto gas = sol->thermo();
    auto kin = sol->kinetics();

    size_t nsp  = gas->nSpecies();

    //--------------- inlet gas state

    //string xin = "H2:1, O2:0.5, N2:1.88";
    string xin = "CH4:1, O2:2, N2:7.52";
    double Tin = 300;
    double P   = 101325;

    vector<double> yin(nsp);
    double         hin;

    gas->setState_TPX(Tin, P, xin);
    gas->getMassFractions(&yin[0]);
    hin = gas->enthalpy_mass();

    //--------------- adiabatic equilibrium gas state

    gas->equilibrate("HP");
    double Tad = gas->temperature();
    vector<double> yad(nsp);
    gas->getMassFractions(&yad[0]);

    cout << endl << "Tad = " << Tad << endl;

    //--------------- PSR object, scaling arrays

    PSR psr(gas, kin, yin, hin, P);

    vector<double> y_tau_scale(nsp+1, 1.0);
    vector<double> f_scale(nsp+1, 1.0);

    //for(int k=0; k<nvar; k++){{{{
    //    ytauscale[k] = ytau[k];
    //    if(ytauscale[k] < 1E-10) ytauscale[k] = 1E-10;
    //    ytauscale[k] = 1.0/ytauscale[k];
    //    //fscale[k] = 1.0/(0.1/tau);
    //}}}}

    //--------------- solve psr

    int    nT   = 500;            // number of T values to solve for
    double Tmin = Tin;
    double Tmax = Tad-0.01;
    vector<double> Tvec(nT);      // temperature values
    for(int i=0; i<nT; i++)
        Tvec[i] = Tmax - (double)(i)/nT * (Tmax - Tmin);

    double taug = 1;           // first guess for tau (large for near equilibrium
    vector<double> y_tau = yad;   // unknown vector: species mass fractions and tau
    y_tau.push_back(taug);

    vector<vector<double>> ys(nT, vector<double>(nsp, 0.0));   // store y vec at each point
    vector<double>         taus(nT);                           // store tau at each point

    ofstream ofile("tau_T.dat");
    ofile << "# tau (s), T (K) " << endl;

    for(int i=0; i<nT; i++) {                                      // loop over each point
        psr.set_T(Tvec[i]);
        psr.solvePSR(y_tau, y_tau_scale, f_scale);                 // solve psr at this point
        ys[i].insert(ys[i].begin(), y_tau.begin(), y_tau.end()-1); // store solution
        taus[i] = y_tau.back();                                    // store solution
        ofile << y_tau.back() << " " << Tvec[i] << endl;
    }

    ofile.close();
    system("ipython plot.py");

    return 0;
}

