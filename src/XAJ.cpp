
#include "math.h"

#include <Rcpp.h>
using namespace Rcpp;

// Run Xinanjiang model
// [[Rcpp::export]]
List XAJrun(NumericVector PREC, NumericVector EVAP, NumericVector parameters,
            NumericVector UH)
{

  /* **************************************************************************
   * Model parameters.
   ************************************************************************* */
  int nrecs = EVAP.length();  // Number of steps

  int i = 0;
  int j = 0;

  double KC = parameters[0];        // 1.  Ratio of potential evap to pan evap
  double IM = parameters[1];        // 2.  Fraction of impermeable area
  double WUM = parameters[2];       // 3.  Soil moisture capacity of upper layer
  double WLM = parameters[3];       // 4.  Soil moisture capacity of lower layer
  double WDM = parameters[4];       // 5.  Soil moisture capacity of deep layer
  double C = parameters[5];         // 6.  Coefficient of deep evap
  double B = parameters[6];         // 7.  Exponent of the soil moisture storage capacity curve

  double SM = parameters[7];        // 8.  Areal mean free water capacity of the surface soil layer
  double EX = parameters[8];        // 9.  Exponent of the free water capacity curve
  double KG = parameters[9];        // 10. outflow coefficients of the free water storage to interflow
  double KI = parameters[10];       // 11. outflow coefficients of the free water storage to groundwater

  double CG = parameters[11];       // 12. recession constant of the lower interflow storage
  double CI = parameters[12];       // 13. recession constant of groundwater storage

  double Area = parameters[13];     // Basin area (km^2)

  double WM = WUM + WLM + WDM;      // Mean water sotrage of the basin
  double WMM = WM * (1. + B) / (1. - IM); // Maximum water storage in the basin

  /** *************************************************************************
   * Output and process variables.
   ************************************************************************* */
  //int nvars = 16;
  //NumericMatrix out(nrecs, nvars);

  NumericVector E_s(nrecs, 0.);
  NumericVector EU_s(nrecs, 0.);
  NumericVector EL_s(nrecs, 0.);
  NumericVector ED_s(nrecs, 0.);
  NumericVector W_s(nrecs, 0.);
  NumericVector WU_s(nrecs, 0.);
  NumericVector WL_s(nrecs, 0.);
  NumericVector WD_s(nrecs, 0.);
  NumericVector R_s(nrecs, 0.);   // Total runoff (mm) of each time step
  NumericVector RS_s(nrecs, 0.);  // Surface runoff (mm) of each time step
  NumericVector RI_s(nrecs, 0.);  // Interflow (mm) of each time step
  NumericVector RG_s(nrecs, 0.);  // Underground runoff (mm) of each time step
  NumericVector Q_s(nrecs, 0.);   // Total runoff (m^3/s) at the outlet of the basin
  NumericVector QS_s(nrecs, 0.);  // Surface runoff (m^3/s) at the outlet of the basin
  NumericVector QI_s(nrecs, 0.);  // Interflow runoff (m^3/s) at the outlet of the basin
  NumericVector QG_s(nrecs, 0.);  // Underground runoff (m^3/s) at the outlet of the basin

  double DeltaT = 24.;    // Hours per time step
  double U = Area/3.6/DeltaT; // Convert runoff from mm to m^3/s

  // Model state variables for each step
  double PE = 0.;  // net prec when > 0; insu evap when < 0
  double Ep = 0.;  // Potential evapotranspiration, KC * EVAP[i]
  double P = 0.;
  double R = 0.;   // Total runoff yield (mm)
  double RB = 0.;  // Runoff yield of the impermeable area
  double RG = 0.;  // Underground runoff yield
  double RI = 0.;  // Interflow runoff yield
  double RS = 0.;  // Surface runoff yield
  double A = 0.;   //
  double E = 0.;   // Total evap
  double EU = 0.;   // Evap at upper soil layer
  double EL = 0.;   // Evap at lower soil layer
  double ED = 0.;   // Evap at deep soil layer

  double S = 0.;
  double FR = 0.;
  double AU = 0.;
  double WU = WUM / 2.;     // Soil moisture (mm) of upper layer
  double WL = WLM * 0.8;    // Soil moisture (mm) of lower layer
  double WD = WDM * 1.0;    // Soil moisture (mm) of deep layer
  double W = WU + WL + WD;  // Soil moisture (mm) of all layer

  double INUL;  // Infiltration from upper layer to lower layer
  double INLD;  // Infiltration from lower layer to deep layer
  double So = 5.0;

  double MS = SM * (1. + EX);
  double FRo = 1 - pow((1. - So / MS), EX);

  /** ************************************************************************
   * Runoff yield
   ************************************************************************* */
  for (i = 0; i < nrecs; i++)
  {
    checkUserInterrupt();

    /**   Evaporation   */
    RB = PREC[i] * IM;     // RB: precipitation of the impermeable area
    P = PREC[i] * (1. - IM);
    Ep = KC * EVAP[i];

    if ((WU + P) >= Ep)
    {
      EU = Ep; EL = 0; ED = 0;
    }
    else if ((WU + P) < Ep)
    {
      EU = WU + P;
      ED = 0.;
      if (WL >= (C * WLM))
      {
        EL = (Ep - EU) * WL / WLM;
      }
      else if (WL >= C * (Ep - EU))
      {
        EL = C * (Ep - EU);
      }
      else if (WL < C * (Ep - EU))
      {
        EL = WL;  ED = C * (Ep - EU) - EL;
      }
      if (ED > WD)
      {
        ED = WD;
      }
    }
    E = EU + EL + ED;
    PE = P - E;

    /**  Infiltration and runoff yeild */
    if (PE <= 0)
    {
      R = 0.00;  W = W + PE;
    }
    else
    {
      A = WMM * (1 - pow((1.0 - W / WM), 1. / (1. + B)));

      // Depth of soil moisture + net prec < maximum soil water storage
      if ((A + PE) < WMM)
      {
        R = PE - (WM - WM * pow((1 - (PE + A) / WMM), (1 + B)) - W);
      }
      // Depth of soil moisture + net prec <= maximum soil water storage
      else
      {
        R = PE - (WM - W);
      }
    }

    // Infiltration
    if (WU + P - EU - R <= WUM)
    {
      WU += P - EU - R;  WL -= EL; WD -= ED;
      INUL = 0.;
    }
    else
    {
      WU = WUM;
      INUL = WU + P - EU - R - WUM;
      if (WL - EL + INUL <= WLM)
      {
        WL = WL - EL + INUL;
        WD = WD - ED;
        INLD = 0.;
      }
      else
      {
        WL = WLM;
        INLD = WL + INUL - EL - WLM;
        if (WD - ED + INLD <= WDM)
          WD = WD - ED + INLD;
        else
          WD = WDM;
      }
    }
    W = WU + WL + WD;
    R += RB;

    // Runoff source division
    if (PE > 0.)
    {
      FR = (R - RB) / PE;
      AU = MS * (1 - pow((1 - So * FRo / FR / SM), 1 / (1 + EX)));
      if (PE + AU < MS)
        RS = FR * (PE + So * FRo / FR - SM + SM * pow((1 - (PE + AU) / MS), 1 + EX));
      else
        RS = FR * (PE + So * FRo / FR - SM);
      S = So * FRo / FR + (R - RS) / FR;
      RI = FR * KI * S;
      RG = FR * KG * S;
      RS += RB;
      R = RS + RI + RG;
      So = S * (1 - KI - KG);
      FRo = FR;
    }
    else
    {
      S = So;
      FR = 1 - pow((1 - S / MS), EX);
      RI = 0.00;
      RG = 0.00;
      So = S * (1 - KI - KG);
      RS = RB;
      R = RS + RI + RG;
      FRo = FR;
    }

    /** Save process variables */
    RS_s[i] = RS;  // Surface runoff (mm) of each time step
    RI_s[i] = RI;  // Interflow (mm) of each time step
    RG_s[i] = RG;  // Underground runoff (mm) of each time step

    E_s[i] = E;
    EU_s[i] = EU;
    EL_s[i] = EL;
    ED_s[i] = ED;
    W_s[i] = W;
    WU_s[i] = WU;
    WL_s[i] = WL;
    WD_s[i] = WD;
    R_s[i] = R;
    RS_s[i] = RS;
    RG_s[i] = RG;
    RI_s[i] = RI;

  }

  /***************************************************************************
   * Routing
   ************************************************************************* */

  // Routing of surface runoff
  int nUH = UH.length();
  int nconv;
  for (i = 0; i < nrecs; i++)
  {
    nconv = nUH;
    if(nUH > i) nconv = i + 1;
    for (j = 0; j < nconv; j++)
    {
      QS_s[i] += RS_s[i-j] * UH[nUH-j-1] * U;
    }
  }

  // Routing of interflow and underground runoff
  for (i = 1; i < nrecs; i++)	{
    QI_s[i] = CI * QI_s[i - 1] + (1.0 - CI) * RI_s[i] * U;
    QG_s[i] = QG_s[i - 1] * CG + RG_s[i] * (1. - CG) * U;
  }

  Q_s = QS_s + QI_s + QG_s;


  // Rcout << "Modeling complete.\n";
  return List::create(E_s, EU_s, EL_s, ED_s, W_s, WU_s, WL_s, WD_s,
                           R_s, RS_s, RG_s, RI_s, Q_s, QS_s, QI_s, QG_s);
}




