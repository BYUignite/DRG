#pragma once

#include <vector>
#include <utility>     // pair, make_pair
#include <memory>      // shared_ptr
#include <set>

#include "cantera/base/Solution.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"

using std::vector;
using std::set;
using std::pair;
using std::shared_ptr;
using Cantera::ThermoPhase;
using Cantera::Kinetics;

////////////////////////////////////////////////////////////////////////////////

class DRG {

    private: 
        shared_ptr<ThermoPhase> gas;
        shared_ptr<Kinetics> kin;
        size_t nsp;
        size_t nrx;
        size_t spPrincipal;    // index of principal species (fuel)
        double eps;

    public:
        set<size_t> spset;     // species set for one state, starting from empty
        set<size_t> spsetU;    // union of species sets over multiple states

    private: 
        Eigen::SparseMatrix<double> RCc;  // reaction coefficients column major; nsp x nrx
        Eigen::SparseMatrix<double,Eigen::RowMajor> RCr;  // reaction coefficients row    major; nsp x nrx

        vector<vector<pair<size_t, vector<size_t>>>> spsprxns; // spsprxn[A][Bloc] is a pair containing a species B and a vector of rxns shared by A and B
        vector<vector<size_t>> sprxns;    // sprxns[A] is a vector of indicies of rxns involving species A
        vector<set<size_t>>    spsp;      // spsp[A] is a set of species B that are paired with A through shared reactions
        vector<set<size_t>>    spSP;      // spSP[A] is the set of species that A depends on (has rAB >= eps)

        //-------------- member functions

        vector<double> get_sumAbsSpRates(const vector<double> &rop);

        double get_rAB(const size_t A, const size_t Bloc, const double sumAbsRatesA, const vector<double> &rop);

        void fill_spset(size_t A);        // recursive function called by DRGspeciesSet

    public:
        void DRGspeciesSet();

        //-------------- constructor functions

        DRG(shared_ptr<ThermoPhase> p_gas, shared_ptr<Kinetics> p_kin, size_t p_spPrincipal, double p_eps);

};
