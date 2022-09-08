#include "psr.h"
#include "streams.h"

#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <utility>
#include <cmath>
#include <set>

using namespace std;
using namespace Cantera;

////////////////////////////////////////////////////////////////////////////////

class DRG {

    public:

    shared_ptr<ThermoPhase> gas;
    shared_ptr<Kinetics> kin;
    size_t nsp;
    size_t nrx;
    size_t spPrincipal;    // index of principal species (fuel)
    double eps;
    set<size_t> spset;

    Eigen::SparseMatrix<double> RCc;  // reaction coefficients column major; nsp x nrx
    Eigen::SparseMatrix<double> RCr;  // reaction coefficients row    major; nsp x nrx

    vector<vector<pair<size_t, vector<size_t>>>> spsprxns; // spsprxn[A][Bloc] is a pair containing a species B and a vector of rxns shared by A and B
    vector<vector<size_t>> sprxns;    // sprxns[A] is a vector of indicies of rxns involving species A
    vector<set<size_t>>    spsp;      // spsp[A] is a set of species B that are paired with A through shared reactions
    vector<set<size_t>>    spSP;      // spSP[A] is the set of species that A depends on (has rAB >= eps)

    //--------------------------------------------------------------------------------

    DRG(shared_ptr<ThermoPhase> p_gas, shared_ptr<Kinetics> p_kin, size_t p_spPrincipal, double p_eps) :
        gas(p_gas), kin(p_kin), spPrincipal(p_spPrincipal), eps(p_eps) {

        nsp = gas->nSpecies();
        nrx = kin->nReactions();

        //--------------------

        RCc = kin->productStoichCoeffs() - kin->reactantStoichCoeffs();
        RCr = Eigen::SparseMatrix<double,Eigen::RowMajor>(RCc);
        RCc.makeCompressed();
        RCr.makeCompressed();

        //-------------------- sprxns

        sprxns.resize(nsp);
        for(size_t A=0; A<nsp; A++) {
            size_t nrxnForSp = RCr.outerIndexPtr()[A+1] - RCr.outerIndexPtr()[A];
            for(size_t i=0; i<nrxnForSp; i++)
                sprxns[A].push_back(RCr.innerIndexPtr()[RCr.outerIndexPtr()[A] + i]);
        }

        //-------------------- spsp

        spsp.resize(nsp);
        for(size_t A=0; A<nsp; A++) {
            size_t nrxnForSp = RCr.outerIndexPtr()[A+1] - RCr.outerIndexPtr()[A];
            for(size_t i=0; i<nrxnForSp; i++) {
                size_t irxn = RCr.innerIndexPtr()[RCr.outerIndexPtr()[A] + i];
                size_t nspInRxn = RCc.outerIndexPtr()[irxn+1] - RCc.outerIndexPtr()[irxn];
                for(size_t k=0; k<nspInRxn; k++) {
                    size_t B = RCc.innerIndexPtr()[RCc.outerIndexPtr()[irxn] + k];
                    spsp[A].insert(B);      // does nothing if B is already in the set
                }
            }
        }

        //-------------------- spsprxns

        spsprxns.resize(nsp);
        for(size_t A=0; A<nsp; A++) {
            for(auto k : spsp[A]) {                    // loop species k paired with A
                spsprxns[A].push_back(make_pair(k, vector<size_t>()));
                for(auto i : sprxns[k]) {              // loop rxns with species k
                    if(RCr.coeff(A, i) != 0)
                        spsprxns[A].back().second.push_back(i);
                }
            }
        }
        

    }

    //--------------------------------------------------------------------------------

    vector<double> get_sumAbsSpRates(const vector<double> &rop){

        //-------- denominator of rAB; for each sp, value = |\nu_{sp,irxn} * rop_irxn|
        //-------- assume gas is already set with the desired chemical state

        vector<double> sumAbsSpRates(nsp, 0.0);
        for(size_t A=0; A<nsp; A++)
            for(auto i : sprxns[A])
                sumAbsSpRates[A] += abs( rop[i] * RCr.coeff(A, i) );
        return sumAbsSpRates;
    }

    //--------------------------------------------------------------------------------

    double get_rAB(const size_t A, const size_t Bloc, const double sumAbsRatesA, const vector<double> &rop) {

        // Bloc is location of species B in spsprxns[A]
        // sumAbsRatesA is deonomiator = sum_i |\nu_A,i * rop_i|

        double rAB = 0.0;
        for(auto i : spsprxns[A][Bloc].second)
            rAB += abs( rop[i] * RCr.coeff(A, i) );

        return rAB / sumAbsRatesA;
    }

    //--------------------------------------------------------------------------------

    void createDRGspeciesSet(){

        vector<double> rop(nrx);
        kin->getNetRatesOfProgress(rop.data());

        vector<double> sumAbsSpRates = get_sumAbsSpRates(rop);

        spSP = spsp;
        for(size_t A=0; A<nsp; A++) {
            for(size_t Bloc=0; Bloc<spsprxns[A].size(); Bloc++) {
                size_t B = spsprxns[A][Bloc].first;
                double rAB = get_rAB(A, Bloc, sumAbsSpRates[A], rop);
                if(rAB < eps) {
                    spSP[A].erase(B);
                }   
            }
        }

        //----------- now turn spSP into a set

        spset.clear();
        fill_spset(spPrincipal);        // recursive
    }

    //--------------------------------------------------------------------------------

    void fill_spset(size_t A) {
        spset.insert(A);
        for(auto k : spSP[A]){
            if(spset.find(k) == spset.end())   // true if k is not in the set
                fill_spset(k);                 // add k and loop over its species
            else
                continue;
        }
    }

};

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

