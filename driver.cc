#include "psr.h"
#include "drg.h"
#include "streams.h"

#include <string>
#include <iostream>

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

    double eps = 0.0;

    DRG drg(gas, kin, gas->speciesIndex("CH4"), eps);

    //--------------- initialize streams

    double P = 101325;
    double T0 = 300;
    double T1 = T0;                  // two streams for composition, but one psr inlet
    string x0 = "O2:2, N2:7.52";
    string x1 = "CH4:1";

    gas->setState_TPX(T0, P, x0);
    double h0 = gas->enthalpy_mass();
    vector<double> y0(nsp);
    gas->getMassFractions(&y0[0]);

    gas->setState_TPX(T1, P, x1);
    double h1 = gas->enthalpy_mass();
    vector<double> y1(nsp);
    gas->getMassFractions(&y1[0]);

    streams strm(gas, P, T0, T1, x0, x1);

    //--------------- storage arrays

    vector<double> yin(nsp);
    double         hin;
    double         Tin;

    vector<double> yad(nsp);
    double         had;
    double         Tad;

    //--------------- PSR object, scaling arrays

    PSR psr(gas, kin);

    vector<double> y_tau_scale(nsp+1, 1.0);
    vector<double> f_scale(nsp+1, 1.0);

    //--------------- solve the psr for each composition

    int nmixf = 10;
    double mixfstart = 0.03; //strm.mixfStoic;
    double mixfend   = 0.07; //mixfstart;
    vector<double> mixfvec(nmixf);
    for(int i=0; i<nmixf; i++)
        mixfvec[i] = mixfstart + (double)(i)/nmixf * (mixfend - mixfstart);

    for(int imixf=0; imixf<nmixf; imixf++) {               // LOOP over each composition

        strm.getMixingState(   mixfvec[imixf], yin, hin, Tin);
        strm.getEquilibrium_HP(mixfvec[imixf], yad, had, Tad);

        psr.setInlet(yin, hin, P);

        cout << endl << imixf << " " << mixfvec[imixf]; cout.flush();

        //--------------- solve psr for each T for given composition

        int    nT   = 200;            // number of T values to solve for; higher values make the solution easier to progress (closer guess)
        double Tmin = Tin;
        double Tmax = Tad - 5.1;      // this can be 0.1 or 0.01 for stoich methane/air, but higher like 5 or more for lean to 0.03 mixf
        vector<double> Tvec(nT);      // temperature values
        for(int i=0; i<nT; i++)
            Tvec[i] = Tmax - (double)(i)/nT * (Tmax - Tmin);

        double taug = 0.01;              // first guess for tau; solver likes a low guess
        vector<double> y_tau = yad;   // unknown vector: species mass fractions and tau
        y_tau.push_back(taug);

        for(int i=0; i<nT; i++) {                          // LOOP over each temperature
            psr.setT(Tvec[i]);
            psr.solvePSR(y_tau, y_tau_scale, f_scale);     // solve psr at this point

            gas->setState_TPY(Tvec[i], P, &y_tau[0]);
            // drg.createDRGspeciesSet();
        }
    }

    return 0;
}

