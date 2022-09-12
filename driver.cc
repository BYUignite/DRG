#include "psr.h"
#include "drg.h"
#include "streams.h"

#include <string>
#include <iostream>

using namespace std;
using namespace Cantera;

////////////////////////////////////////////////////////////////////////////////

int main() {

    //--------------- user inputs

    string mechName   = "gri30.yaml";
    string skMechName = "skel.yaml";
    double eps = 0.3;
    vector<string> spPrincipal = {"CH4"};
    vector<string> spExtra     = {"N2"};

    double P = 101325;
    double T0 = 300;
    string x0 = "O2:2, N2:7.52";
    string x1 = "CH4:1";

    int nmixf = 11;
    double mixfstart = 0.03; //strm.mixfStoic;
    double mixfend   = 0.07; //mixfstart;

    int    nT   = 201;            // number of T values to solve for; higher values make the solution easier to progress (closer guess)
    double Tmin = 1000;
    double TmaxDelta = -5.1;      // this can be 0.1 or 0.01 for stoich methane/air, but higher like 5 or more for lean to 0.03 mixf

    double taug = 0.01;           // first guess for tau; solver likes a low guess

    //--------------- initialize cantera

    double T1 = T0;                  // two streams for composition, but one psr inlet

    auto sol = newSolution(mechName, "", "None");
    auto gas = sol->thermo();
    auto kin = sol->kinetics();

    size_t nsp  = gas->nSpecies();

    DRG drg(gas, kin, spPrincipal, spExtra, eps);

    //--------------- initialize streams

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

    vector<double> yad(nsp);
    double         had;
    double         Tad;

    //--------------- PSR object, scaling arrays

    PSR psr(gas, kin);

    vector<double> y_tau_scale(nsp+1, 1.0);
    vector<double> f_scale(nsp+1, 1.0);

    //--------------- solve the psr for each composition

    vector<double> mixfvec(nmixf);
    for(int i=0; i<nmixf; i++)
        mixfvec[i] = mixfstart + (double)(i)/(nmixf-1) * (mixfend - mixfstart);

    cout << endl << "Solving full PSR S-curve for the following mixture fractions: ";

    double Tdmb;

    for(int imixf=0; imixf<nmixf; imixf++) {               // LOOP over each composition

        strm.getMixingState(   mixfvec[imixf], yin, hin, Tdmb);
        strm.getEquilibrium_HP(mixfvec[imixf], yad, had, Tad);

        psr.setInlet(yin, hin, P);

        cout << endl << mixfvec[imixf];

        //--------------- solve psr for each T for given composition

        double Tmax = Tad + TmaxDelta;      // this can be 0.1 or 0.01 for stoich methane/air, but higher like 5 or more for lean to 0.03 mixf
        vector<double> Tvec(nT);      // temperature values
        for(int i=0; i<nT; i++)
            Tvec[i] = Tmax - (double)(i)/(nT-1) * (Tmax - Tmin);

        vector<double> y_tau = yad;   // unknown vector: species mass fractions and tau
        y_tau.push_back(taug);

        for(int i=0; i<nT; i++) {                          // LOOP over each temperature
            psr.setT(Tvec[i]);
            psr.solvePSR(y_tau, y_tau_scale, f_scale);     // solve psr at this point

            gas->setState_TPY(Tvec[i], P, &y_tau[0]);
            drg.DRGspeciesSet();
        }
    }

    //--------------- output the skeletal species

    cout << endl;
    cout << endl << "Cantera mechanism name: " << sol->name();
    cout << endl << "DRG tolerance: " << eps;

    //--------------- create skeletal mechanism

    drg.writeSkeletalMechanism(mechName, skMechName);

    return 0;
}

