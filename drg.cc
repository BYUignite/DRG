#include "drg.h"

#include <cmath>
#include <iostream>

using namespace std;
using namespace Cantera;

////////////////////////////////////////////////////////////////////////////////

DRG::DRG(shared_ptr<ThermoPhase> p_gas, shared_ptr<Kinetics> p_kin, 
         size_t p_spPrincipal, double p_eps) :
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

////////////////////////////////////////////////////////////////////////////////

vector<double> DRG::get_sumAbsSpRates(const vector<double> &rop){

    //-------- denominator of rAB; for each sp, value = |\nu_{sp,irxn} * rop_irxn|
    //-------- assume gas is already set with the desired chemical state

    vector<double> sumAbsSpRates(nsp, 0.0);
    for(size_t A=0; A<nsp; A++)
        for(auto i : sprxns[A])
            sumAbsSpRates[A] += abs( rop[i] * RCr.coeff(A, i) );
    return sumAbsSpRates;
}

////////////////////////////////////////////////////////////////////////////////

double DRG::get_rAB(const size_t A, const size_t Bloc, const double sumAbsRatesA, const vector<double> &rop) {

    // Bloc is location of species B in spsprxns[A]
    // sumAbsRatesA is deonomiator = sum_i |\nu_A,i * rop_i|

    double rAB = 0.0;
    for(auto i : spsprxns[A][Bloc].second)
        rAB += abs( rop[i] * RCr.coeff(A, i) );

    return rAB / sumAbsRatesA;
}

////////////////////////////////////////////////////////////////////////////////

void DRG::createDRGspeciesSet(){

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

////////////////////////////////////////////////////////////////////////////////

void DRG::fill_spset(size_t A) {
    spset.insert(A);
    for(auto k : spSP[A]){
        if(spset.find(k) == spset.end())   // true if k is not in the set
            fill_spset(k);                 // add k and loop over its species
        else
            continue;
    }
}

