#include "drg.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <regex>

#include <cstdlib>

using namespace std;
using namespace Cantera;

////////////////////////////////////////////////////////////////////////////////

DRG::DRG(shared_ptr<ThermoPhase> p_gas, shared_ptr<Kinetics> p_kin, 
         vector<string> p_spPrincipal, vector<string> p_Ex, double p_eps) :
    gas(p_gas), kin(p_kin), spPrincipal(p_spPrincipal), spExtra(p_Ex), eps(p_eps) {

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
                if(B != A)
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

void DRG::DRGspeciesSet(){

    vector<double> rop(nrx);
    kin->getNetRatesOfProgress(rop.data());

    vector<double> sumAbsSpRates = get_sumAbsSpRates(rop);

    spSP = spsp;
    for(size_t A=0; A<nsp; A++) {
        for(size_t Bloc=0; Bloc<spsprxns[A].size(); Bloc++) {
            size_t B = spsprxns[A][Bloc].first;
            double rAB = get_rAB(A, Bloc, sumAbsSpRates[A], rop);
            //cout << endl << gas->speciesName(A) << " " << gas->speciesName(B) << " " << rAB;
            if(rAB < eps) {
                spSP[A].erase(B);
            }   
        }
    }

    //----------- now turn spSP into a set

    spset.clear();
    for(auto spName : spPrincipal) {
        size_t A = gas->speciesIndex(spName);
        if(spset.find(A) == spset.end())
            fill_spset(A);                   // recursive
    }

    //----------- merge spset with: spsetU from previous runs, = union

    spsetU.merge(spset);                             // c++17 feature
    //set_union(spset.begin(), spset.end(),          // or this nastiness, but sorts too
    //          spsetU.begin(), spsetU.end(), 
    //          inserter(spsetU, spsetU.begin()));

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

////////////////////////////////////////////////////////////////////////////////

void DRG::DRGreactionsSet(){

    for(auto A : spsetU) {                           // loop over skeletal species A
        for(auto irxn : sprxns[A]) {                 // loop over reactions that include A
            if(rxsetU.find(irxn) != rxsetU.end())    // don't process reaction if it's already included
                continue;
            bool keepRxn = true;
            size_t nspInRxn = RCc.outerIndexPtr()[irxn+1] - RCc.outerIndexPtr()[irxn];
            for(size_t kk=0; kk<nspInRxn; kk++) {    // loop over species k in reaction irxn
                size_t k = RCc.innerIndexPtr()[RCc.outerIndexPtr()[irxn] + kk];
                if(k==A) continue;
                if(spsetU.find(k) == spsetU.end()) { // rxn has a species not in the skeletal list, don't include rxn
                    keepRxn = false;
                    break;
                }
            }
            if(keepRxn)
                rxsetU.insert(irxn);
        }
    }

    // ISSUE: reaction with species on both sides of an equation is not reported in the coefficient matrices
    // Use regular expression to find offenders and strip out the relevant reaction.
    // This might result in some inconsistency with the mechanism and the tolerance chosen,
    // but that tolerance is not particularly important; the skeletal mechanism will require validation anyway.

    set<size_t> nonSkelSpecies;           // species in the detailed mechanism that are not in the skeletal mechanism
    for(size_t k=0; k<nsp; k++)
        if(spsetU.find(k) == spsetU.end())
            nonSkelSpecies.insert(k);
    
    string rxn;
    string sp;
    regex rs, re, rm;
    vector<size_t> toerase;
    for(auto k : nonSkelSpecies) {
        sp = gas->speciesName(k);
        rs = ("^" + sp + " ");
        re = (" " + sp + "$");
        rm = (" " + sp + " ");
        for(auto i : rxsetU) {
            rxn = kin->reactionString(i);
            if(regex_search(rxn, rs) || regex_search(rxn, re) || regex_search(rxn, rm)) {
                toerase.push_back(i);
            }
        }
    }
    for(auto k : toerase)
        rxsetU.erase(k);
}

////////////////////////////////////////////////////////////////////////////////
// fdetMech is the name of the detailed mechanism that everything is based on (*.yaml)
// fsklMech is the name of the new skeletal mechanism file (*.yaml)

void DRG::writeSkeletalMechanism(string fdetMech, string fsklMech) {

    //----------- add extra species to spsetU

    for(auto spName : spExtra)
        spsetU.insert(gas->speciesIndex(spName));

    //----------- display species

    cout << endl << "# species: " << spsetU.size();
    for(auto k : spsetU)
        cout << endl << "    " << gas->speciesName(k);
    cout << endl;

    //----------- get the reaction set corresponding to spsetU

    DRGreactionsSet();

    cout << endl << "# reactions: " << rxsetU.size() << endl;

    //-----------

    auto rootNd  = AnyMap::fromYamlFile(fdetMech);
    auto phaseNd = rootNd["phases"].getMapWhere("name", "");

    //----------- reset species array in the phase node

    auto spVec   = phaseNd["species"].as<std::vector<std::string>>();
    for(int k=nsp-1; k>=0; k--){               // count down so deletions dont screw up the iteration
        if(spsetU.find(k) == spsetU.end())     // species k not present in skeletal mechanism
            spVec.erase(spVec.begin()+k);
    }
    rootNd["phases"].getMapWhere("name","")["species"] = spVec;   // udate node with skel species

    //----------- reset species thermo and transport information

    auto spInfoVec = rootNd["species"].asVector<AnyMap>();
    for(int k=spInfoVec.size()-1; k>=0; k--){  // count down as above; loop spInfoVec not orig. gas species list, since latter may be smaller
        string spname = spInfoVec[k]["name"].asString();
        size_t isp = gas->speciesIndex(spname);
        if(spsetU.find(isp) == spsetU.end())   // species k not present in skeletal mechanism
            spInfoVec.erase(spInfoVec.begin() + k);
    }
    rootNd["species"] = spInfoVec;             // update node with skel species info

    //----------- reset reactions

    auto rxInfoVec = rootNd["reactions"].asVector<AnyMap>();
    for(int i=rxInfoVec.size()-1; i>=0; i--){  // count down as above
        if(rxsetU.find(i) == rxsetU.end())     // reaction i not present in skeletal mechanism
            rxInfoVec.erase(rxInfoVec.begin()+i);
    }
    rootNd["reactions"] = rxInfoVec;           // update node with skel reactions info

    //----------- third body efficiencies

    rxInfoVec = rootNd["reactions"].asVector<AnyMap>();
    for(size_t irxn=0; irxn<rxInfoVec.size(); irxn++){
        if(rxInfoVec[irxn].hasKey("efficiencies")){
            auto effs = rxInfoVec[irxn]["efficiencies"].asMap<double>();
            vector<string> toerase(0);
            for(auto eff = effs.begin(); eff != effs.end(); eff++) {
                string spName = eff->first;
                if(spsetU.find(gas->speciesIndex(spName)) == spsetU.end())
                    toerase.push_back(spName);
            }
            for(auto spName : toerase)
                effs.erase(spName);
            rxInfoVec[irxn]["efficiencies"] = effs;
        }
    }
    rootNd["reactions"] = rxInfoVec;

    //----------- write the mechanism

    string mechString = rootNd.toYamlString();
    ofstream ofile(fsklMech);
    ofile << mechString;
    ofile.close();

}
