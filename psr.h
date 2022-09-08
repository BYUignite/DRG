#include <kinsol/kinsol.h>             /* access to KINSOL func., consts. */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector       */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix       */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype */

#include "cantera/base/Solution.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"

#include <memory>
#include <vector>

using std::vector;
using std::shared_ptr;
using namespace Cantera;

////////////////////////////////////////////////////////////////////////////////

static int func(N_Vector yT, N_Vector f, void *user_data);

////////////////////////////////////////////////////////////////////////////////

class PSR {

    public:

        shared_ptr<ThermoPhase> gas;
        shared_ptr<Kinetics>    kin;

        vector<double>          yin;
        double                  hin;
        double                  P;
        double                  T;

        size_t                  neq;

        //-----------------

        PSR(shared_ptr<ThermoPhase> p_gas, 
            shared_ptr<Kinetics>    p_kin)   : 
                gas(p_gas), kin(p_kin) {
            neq = gas->nSpecies() + 1; 
        }

        PSR(shared_ptr<ThermoPhase> p_gas, 
            shared_ptr<Kinetics>    p_kin,
            vector<double>         &p_yin,
            double                  p_hin,
            double                  p_P)     : 
                gas(p_gas), kin(p_kin),
                yin(p_yin), hin(p_hin), P(p_P) {
            neq = gas->nSpecies() + 1; 
        }

        void setInlet(vector<double> &p_yin, double p_hin, double p_P) { 
            yin = p_yin;
            hin = p_hin;
            P   = p_P;
        }

        void setT(double p_T) { T=p_T; }

        //-----------------

        void solvePSR(vector<double> &ytau, 
                      vector<double> &ytauscale, vector<double> &fscale, 
                      double ftol=1.0E-5, double stol=1.0E-5) {

            // ytau: guess on input; solution on output

            SUNContext sun;
            int rv = SUNContext_Create(NULL, &sun);

            N_Vector ytausun  = N_VNew_Serial(neq, sun);  // solution
            N_Vector scl_ytau = N_VNew_Serial(neq, sun);  // solution scaling array
            N_Vector scl_f    = N_VNew_Serial(neq, sun);  // function scaling array

            N_Vector cstrt    = N_VNew_Serial(neq, sun);  // function scaling array
            
            for(size_t k=0; k<neq; k++) {
                NV_Ith_S(ytausun, k)  = ytau[k];
                NV_Ith_S(scl_ytau, k) = ytauscale[k];
                NV_Ith_S(scl_f,  k)   = fscale[k];
                NV_Ith_S(cstrt,  k)   = 0.0;
            }
            NV_Ith_S(cstrt,  neq-1)   = 1.0;

            void * kmem = KINCreate(sun);

            rv = KINInit(kmem, func, ytausun);
            rv = KINSetUserData(kmem, this);
            rv = KINSetConstraints(kmem, cstrt);
            rv = KINSetScaledStepTol(kmem, stol);
            rv = KINSetFuncNormTol(  kmem, ftol);

            SUNMatrix       J  = SUNDenseMatrix(neq, neq, sun);   // linear solver matrix J
            SUNLinearSolver LS = SUNLinSol_Dense(ytausun, J, sun);  // set linear solver
            rv = KINSetLinearSolver(kmem, LS, J);                 // associate matrix J with solver LS

            //---------------------------

            int lnsearch = KIN_LINESEARCH;                  // KIN_NONE, or KIN_LINESEARCH
            int exNewt_modNewt  = 1;                        // =1 Exact Newton; =0 Modified Newton
            rv = KINSetMaxSetupCalls(kmem, exNewt_modNewt);
            rv = KINSol(kmem, ytausun, lnsearch, scl_ytau, scl_f);

            //---------------------------

            for(size_t k=0; k<neq; k++)
                ytau[k] = NV_Ith_S(ytausun,k);

            //---------------- Free memory

            N_VDestroy(ytausun);
            N_VDestroy(scl_ytau);
            N_VDestroy(scl_f);
            N_VDestroy(cstrt);
            KINFree(&kmem);
            SUNLinSolFree(LS);
            SUNMatDestroy(J);
            SUNContext_Free(&sun);
        }
};

////////////////////////////////////////////////////////////////////////////////

static int func(N_Vector ytauvec, N_Vector fvec, void *user_data) {

    PSR *psr = static_cast<PSR *>(user_data);

    realtype *ytau = N_VGetArrayPointer(ytauvec);
    realtype *f  = N_VGetArrayPointer(fvec);

    //psr->gas->setMassFractions(ytau);
    psr->gas->setMassFractions_NoNorm(ytau);
    psr->gas->setState_TP(psr->T, psr->P);
    double rho = psr->gas->density();
    vector<double> rr(psr->neq-1);
    psr->kin->getNetProductionRates(&rr[0]);

    for(size_t k=0; k<psr->neq-1; k++) {
        f[k] = (psr->yin[k] - ytau[k])/ytau[psr->neq-1] + rr[k]*psr->gas->molecularWeight(k)/rho;
    }
    f[psr->neq-1] = psr->gas->enthalpy_mass() - psr->hin;

    return 0;
}
