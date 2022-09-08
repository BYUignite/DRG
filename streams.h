#pragma once

#include "cantera/base/Solution.h"
#include "cantera/thermo.h"

#include <memory>
#include <vector>
#include <string>

using std::vector;
using std::shared_ptr;
using std::string;
using namespace Cantera;

////////////////////////////////////////////////////////////////////////////////

class streams {

    public:

    //////////////////// DATA MEMBERS //////////////////////

        shared_ptr<ThermoPhase> gas;

        double         P;               ///< pressure (Pa)
        double         T0;              ///< stream mixf=0 temperature (K)
        double         T1;              ///< stream mixf=1 temperature (K)
        double         h0;              ///< stream mixf=0 enthalpy (J/kg)
        double         h1;              ///< stream mixf=1 enthalpy (J/kg)
        vector<double> y0;              ///< stream mixf=0 composition vector
        vector<double> y1;              ///< stream mixf=1 composition vector

        double         mixfStoic;       ///< stoichiometric mixture fraction

        int            nspc;            ///< number of species in gas mechanism

                                        /// \fun{\text{mixf} = \frac{\beta-\beta_0}{\beta_1-\beta_0}}
        double         beta0;           ///< mixf = (beta-beta0) / (beta1-beta0)
                                        ///< \fun{\text{mixf} = \frac{\beta-\beta_0}{\beta_1-\beta_0}}
        double         beta1;           ///< mixf = (beta-beta0) / (beta1-beta0)

        vector<double> gCHON;           ///< gammas, as in beta = sum_i (y_i*gamma_i)

    //////////////////// MEMBER FUNCTIONS /////////////////

        void getProdOfCompleteComb(const double mixf,
                                   vector<double> &ypcc,
                                   double &hpcc,
                                   double &Tpcc);

        void getEquilibrium_HP(const double mixf,
                                     vector<double> &yeq,
                                     double &heq,
                                     double &Teq);

        void getEquilibrium_TP(const double mixf,
                                     double Teq,
                                     vector<double> &yeq,
                                     double &heq);

        void getMixingState(const double mixf,
                            vector<double> &ymix,
                            double &hmix,
                            double &Tmix);

        void setGasMixingState(const double mixf);

        double getMixtureFraction(const double *y,
                                  const bool doBeta01=false);

    private:

        void setStoicMixf();
        vector<double> setElementMassFracs(const double *y);
        vector<double> setElementMoleFracs(const double *y);
        vector<double> getElementMoles(const double *x,
                                       double &nOnotFromO2,
                                       double &nHnotFromH2O,
                                       double &nCnotFromCO2);

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        streams(shared_ptr<ThermoPhase> p_gas, 
                double p_P, double p_T0, double p_T1,
                string x0, string x1);

        ~streams(){}
};
