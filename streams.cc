
#include "streams.h"
#include <cassert>
#include <iostream>

using namespace std;

///////////////////////////////////////////////////////////////////////////////
/** Constructor
 *
 * @param p_gas \input pointer to cantera ThermoPhase object (gas)
 * @param p_P   \input system pressure (Pa)
 * @param p_T0  \input temperature of stream 0
 * @param p_T1  \input temperature of stream 1
 * @param x0    \input string mole or mole fraction composition
 * @param x1    \input string mole or mole fraction composition
 */
streams::streams(shared_ptr<ThermoPhase> p_gas, 
                 double p_P, double p_T0, double p_T1, 
                 string x0, string x1) : 
    gas(p_gas), P(p_P), T0(p_T0), T1(p_T1) {

    nspc = gas->nSpecies();

    y0   = vector<double>(nspc, 0.0);
    y1   = vector<double>(nspc, 0.0);

    gas->setState_TPX(T0, P, x0);
    gas->getMassFractions(&y0[0]);
    h0     = gas->enthalpy_mass();

    gas->setState_TPX(T1, P, x1);
    gas->getMassFractions(&y1[0]);
    h1     = gas->enthalpy_mass();

    gCHON.resize(4);
    gCHON[0] =  2.0/gas->atomicWeight(gas->elementIndex("C"));
    gCHON[1] =  0.5/gas->atomicWeight(gas->elementIndex("H"));
    gCHON[2] = -1.0/gas->atomicWeight(gas->elementIndex("O"));
    gCHON[3] =  0.0;

    //---------------- stoich mixf

    setStoicMixf();

    getMixtureFraction(&y0[0], true);    // set beta0 and beta1

}
///////////////////////////////////////////////////////////////////////////////
/** Computes the temperature, enthalpy, and composition of mixing among streams.
 *
 *  @param mixf \input mixture fraction, defines elemental composition.
 *  @param ymix \output mass fractions of products of complete combustion.
 *  @param hmix \output enthalpy of products of complete combustion.
 *  @param Tmix \output temperature of products of complete combustion.
 */
void streams::getMixingState(const double mixf, vector<double> &ymix,
                             double &hmix, double &Tmix){

    hmix = h1*mixf + h0*(1.0-mixf);
    for(int k=0; k<nspc; k++)
        ymix[k] = y1[k]*mixf + y0[k]*(1.0-mixf);

    gas->setMassFractions( &ymix[0] );
    gas->setState_HP(hmix, P, 1.E-10);    // get temperature as Tadiabatic

    Tmix = gas->temperature();
}

///////////////////////////////////////////////////////////////////////////////
/** Set the gas state based on the mixture fraction
 *
 *  @param mixf \input mixture fraction, defines elemental composition.
 */
void streams::setGasMixingState(const double mixf){

    double hmix = h1*mixf + h0*(1.0-mixf);
    vector<double> ymix(nspc);
    for(int k=0; k<nspc; k++)
        ymix[k] = y1[k]*mixf + y0[k]*(1.0-mixf);

    gas->setMassFractions( ymix.data() );
    gas->setState_HP(hmix, P, 1.E-10);
}

///////////////////////////////////////////////////////////////////////////////
/** Computes the temperature, enthalpy, and composition of equilibrium at the given mixf.
 *
 *  @param mixf \input mixture fraction, defines elemental composition.
 *  @param yeq \output mass fractions
 *  @param heq \output enthalpy
 *  @param Teq \output temperature
 */
void streams::getEquilibrium_HP(const double mixf, vector<double> &yeq,
                                      double &heq, double &Teq){

    //---------- Compute the mixing mass fractions and enthalpy

    yeq.resize(nspc);
    for(int k=0; k<nspc; k++)
        yeq[k] = y1[k]*mixf + y0[k]*(1.0-mixf);
    heq = h1*mixf + h0*(1.0-mixf);

    gas->setState_PY(P, &yeq[0]);
    gas->setState_HP(heq, P, 1.E-10);

    gas->equilibrate("HP");
    gas->getMassFractions(&yeq[0]);

    Teq = gas->temperature();

}

///////////////////////////////////////////////////////////////////////////////
/** Computes the enthalpy, and composition of equilibrium at the given mixf for given T
 *
 *  @param mixf \input mixture fraction, defines elemental composition.
 *  @param Teq \input temperature
 *  @param yeq \output mass fractions
 *  @param heq \output enthalpy
 */
void streams::getEquilibrium_TP(const double mixf, double Teq, 
                                      vector<double> &yeq, double &heq ){

    //---------- Compute the mixing mass fractions and enthalpy

    yeq.resize(nspc);
    for(int k=0; k<nspc; k++)
        yeq[k] = y1[k]*mixf + y0[k]*(1.0-mixf);

    gas->setState_TPY(Teq, P, &yeq[0]);

    gas->equilibrate("TP");
    gas->getMassFractions(&yeq[0]);

    heq = gas->enthalpy_mass();

}

///////////////////////////////////////////////////////////////////////////////
/** Computes the temperature, enthalpy, and composition of complete combustion at the given mixf.
 *  For nonpremixed flames (don't do anything funny, like have oxygen in the fuel stream)
 *
 *  @param mixf \input mixture fraction, defines elemental composition.
 *  @param ypcc \output mass fractions of products of complete combustion.
 *  @param hpcc \output enthalpy of products of complete combustion.
 *  @param Tpcc \output temperature of products of complete combustion.
 */
void streams::getProdOfCompleteComb(const double mixf, vector<double> &ypcc,
                                    double &hpcc, double &Tpcc){

    //---------- Compute the mixing mass fractions and enthalpy

    ypcc.resize(nspc);
    double d1 = 1.0-mixf;
    for(int k=0; k<nspc; k++) {
        ypcc[k] = y1[k]*mixf + y0[k]*d1;
    }
    hpcc = h1*mixf + h0*d1;

    //--------- Set gas and element indicicies

    int iC = gas->elementIndex("C");
    int iH = gas->elementIndex("H");
    //int iO = gas->elementIndex("O"); // !!!!!  currently unusedvariable
    int iN = gas->elementIndex("N");
    int iCO2 = gas->speciesIndex("CO2");
    int iH2O = gas->speciesIndex("H2O");
    int iN2  = gas->speciesIndex("N2");
    int iO2  = gas->speciesIndex("O2");

    //---------- Set ypcc as the mixing mole fractions: Take a basis of one mole:
    // now we are working in moles
    // elemM are moles of each element
    // when stoic:
    // CxHyNz   + (x+y/4)O2  ==>  (z/2)N2 + (x)CO2 + (y/2)H2O
    // otherwise:
    // CxHyNz   + (beta)O2   ==>  (z/2)N2 + (x)CO2 + (y/2)H2O
    // Note this allows fuels with nitrogen and oxygen

    gas->setMassFractions( &ypcc[0] );
    gas->getMoleFractions( &ypcc[0] );

    double nOnotFromO2  = 0.0;
    double nHnotFromH2O = 0.0;
    double nCnotFromCO2 = 0.0;
    vector<double> elemM = getElementMoles( &ypcc[0], nOnotFromO2,
            nHnotFromH2O, nCnotFromCO2 );

    double x    = elemM[iC];
    double y    = elemM[iH];
    double z    = elemM[iN];
    double beta = elemM[gas->elementIndex("O")] * 0.5;        // moles of O as O2

    //--------------------------------------------------------------------------

    if(mixf < mixfStoic) {                        // lean: burn all fuel, leftover O2

        for(int k=0; k<nspc; k++)
            ypcc[k] = 0.0;

        if(iCO2 > 0)                              // CO2 is not in the H2 mechs
            ypcc[iCO2] = x;
        ypcc[iH2O] = y*0.5;
        ypcc[iN2]  = z*0.5;
        ypcc[iO2]  = beta - (x+y/4.0);

    }
    else{                                         // rich: burn all O2, leftover fuel

        //double eta = beta/(x+y/4.0); // extent of reaction
        double eta = (beta-nOnotFromO2/2)/(x+y/4-nOnotFromO2/2); // extent of reaction
        if(eta > 1.0)
            cout << endl << "eta > 1.0" << endl;
        d1 = 1.0-eta;                            // fraction of fuel unburnt

        for(int k=0; k<nspc; k++)
            ypcc[k] *= d1;                       // initialize everything then correct
        if(iCO2 > 0)                             // CO2 is not in the H2 mechs
            ypcc[iCO2] = (x-nCnotFromCO2) + nCnotFromCO2*eta;       // init + formed
        ypcc[iH2O] = (y-nHnotFromH2O)*0.5 + nHnotFromH2O*0.5*eta;   // init + formed
        ypcc[iN2]  = z*0.5;
        ypcc[iO2]  = 0.0;

    }

    //--------------------------------------------------------------------------

    double sum = 0.0;                       // normalize moles
    for(int k=0; k<nspc; k++)
        sum += ypcc[k];
    assert(sum != 0.0);
    for(int k=0; k<nspc; k++)
        ypcc[k] /= sum;
    gas->setMoleFractions( &ypcc[0] );      // set mole fractions
    gas->getMassFractions( &ypcc[0] );      // set ypcc as mass fractions


    //--------------------------------------------------------------------------

    gas->setState_HP(hpcc, P, 1.E-10);    // get temperature as Tadiabatic
    //gas->setState_HP(hpcc, P);    // get temperature as Tadiabatic
    Tpcc = gas->temperature();

}

///////////////////////////////////////////////////////////////////////////////
/** Set the stoichiometric mixture fraction using Bilger's definition */

void streams::setStoicMixf() {

    vector<double> x_c(nspc,0.0);

    double mc0, mc1, mo0, mo1, mh0, mh1;

    vector<double> elemMassFrac0 = setElementMassFracs(&y0[0]);
    vector<double> elemMassFrac1 = setElementMassFracs(&y1[0]);

    mc0 = elemMassFrac0[gas->elementIndex("C")]/
        gas->atomicWeight(gas->elementIndex("C"));
    mc1 = elemMassFrac1[gas->elementIndex("C")]/
        gas->atomicWeight(gas->elementIndex("C"));
    mh0 = elemMassFrac0[gas->elementIndex("H")]/
        gas->atomicWeight(gas->elementIndex("H"));
    mh1 = elemMassFrac1[gas->elementIndex("H")]/
        gas->atomicWeight(gas->elementIndex("H"));
    mo0 = elemMassFrac0[gas->elementIndex("O")]/
        gas->atomicWeight(gas->elementIndex("O"));
    mo1 = elemMassFrac1[gas->elementIndex("O")]/
        gas->atomicWeight(gas->elementIndex("O"));

    mixfStoic = (2.0*mc0 + 0.5*mh0 - mo0) /
        (mo1-mo0 + 2.0*(mc0-mc1) + 0.5*(mh0-mh1));

}

///////////////////////////////////////////////////////////////////////////////
/** Sets the elements to have the correct Mass Fractions based on the specified array.
 *  @param y \input mass fraction array to use to get corresponding element fractions.
 *  @return vector of element mass fractions.
 */

vector<double> streams::setElementMassFracs(const double *y) {


    vector<double> atomArr(gas->nElements());
    vector<double> elemMassFrac(gas->nElements(), 0.0);
    double sum = 0.0;

    gas->setMassFractions( &y[0] );

    for(int k=0; k<nspc; k++) {
        sum=0.0;
        gas->getAtoms(k, &atomArr[0]);           // [nelements] in sp k
        for(int m=0; m<(int)atomArr.size(); m++)
            sum += atomArr[m] * gas->atomicWeight(m);
        for(int m=0; m<(int)atomArr.size(); m++)
            elemMassFrac[m] += y[k] * atomArr[m]/sum * gas->atomicWeight(m);
                              // is * mass frac of elem in sp
    }

    return elemMassFrac;

}

///////////////////////////////////////////////////////////////////////////////
/** Sets the elements to have the correct Mole Fractions based on the specified array.
 *  @param y \input mass fraction array to use to get corresponding element fractions.
 *  @return vector of element mole fractions.
 */

vector<double> streams::setElementMoleFracs(const double *y) {


    vector<double> atomArr(gas->nElements());
    vector<double> elemMoleFrac(gas->nElements(), 0.0);

    gas->setMassFractions( &y[0] );

    for(int k=0; k<nspc; k++) {
        gas->getAtoms(k, &atomArr[0]);           // [nelements] in sp k
        for(int m=0; m<(int)atomArr.size(); m++)
            elemMoleFrac[m] += gas->moleFraction(k) * atomArr[m];
    }
    double sum = 0.0;
    for(int m=0; m<(int)atomArr.size(); m++)
        sum += elemMoleFrac[m];
    assert(sum != 0.0);
    for(int m=0; m<(int)atomArr.size(); m++)
        elemMoleFrac[m] /= sum;

    return elemMoleFrac;

}

///////////////////////////////////////////////////////////////////////////////
/** Get amount of moles for each element.
 *  @param x \input pointer to vector of species mole fractions.
 *  @param nOnotFromO2 \input number of moles of oxygen not from O2 (oxygen in the base fuel).
 *  @param nHnotFromH2O \input number of moles of hydrogen not from H2O.
 *  @param nCnotFromCO2 \input number of moles of carbon not from CO2.
 *  @return vector of element moles.
 */

vector<double> streams::getElementMoles(const double *x,
                                        double &nOnotFromO2,
                                        double &nHnotFromH2O,
                                        double &nCnotFromCO2) {


    vector<double> atomArr(gas->nElements());
    vector<double> elemM(gas->nElements(), 0.0);
    int iO2  = gas->speciesIndex("O2");
    int iO   = gas->elementIndex("O");
    int iCO2 = gas->speciesIndex("CO2");
    int iC   = gas->elementIndex("C");
    int iH2O = gas->speciesIndex("H2O");
    int iH   = gas->elementIndex("H");

    for(int k=0; k<nspc; k++) {
        gas->getAtoms(k, &atomArr[0]);           // [nelements] in sp k
        for(int m=0; m<(int)atomArr.size(); m++)
            elemM[m] += x[k] * atomArr[m];
        if(k != iO2)  nOnotFromO2  += atomArr[iO] * x[k];
        if(k != iCO2) nCnotFromCO2 += atomArr[iC] * x[k];
        if(k != iH2O) nHnotFromH2O += atomArr[iH] * x[k];
    }
    return elemM;

}

///////////////////////////////////////////////////////////////////////////////
/**Compute the mixture fraction from the mass fractions using Bilger's mixf.
 * Set doBeta01=true on first call to initialize members beta0, beta1.
 * Later calls of this function only use the first parameter.
 *
 * @param y \input vector of species mass fractions.
 * @param doBeta01 \input flag=true on first call to set members beta0, beta1.
 * @return mixture fraction
 */

double streams::getMixtureFraction(const double *y, const bool doBeta01) {

    vector<double> elemMF;
    double         beta;

    elemMF = setElementMassFracs(y);
    beta   = gCHON[0] * elemMF[gas->elementIndex("C")] +
             gCHON[1] * elemMF[gas->elementIndex("H")] +
             gCHON[2] * elemMF[gas->elementIndex("O")] +
             gCHON[3] * elemMF[gas->elementIndex("N")];

    if(doBeta01) {
        elemMF = setElementMassFracs(&y0[0]);
        beta0  = gCHON[0] * elemMF[gas->elementIndex("C")] +
                 gCHON[1] * elemMF[gas->elementIndex("H")] +
                 gCHON[2] * elemMF[gas->elementIndex("O")] +
                 gCHON[3] * elemMF[gas->elementIndex("N")];

        elemMF = setElementMassFracs(&y1[0]);
        beta1  = gCHON[0] * elemMF[gas->elementIndex("C")] +
                 gCHON[1] * elemMF[gas->elementIndex("H")] +
                 gCHON[2] * elemMF[gas->elementIndex("O")] +
                 gCHON[3] * elemMF[gas->elementIndex("N")];
    }

    if( beta1-beta0 == 0 )
        return 0.0;            // to avoid division by zero
    else
        return (beta - beta0)/(beta1 - beta0);

}

