#include "KLNModel.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "Regge96.h"

#ifdef ALPHA
#define abs  fabs
inline double max(double a, double b) { return a>=b ? a:b;}
inline double min(double a, double b) { return a<=b ? a:b;}
#endif

using namespace std;

const double hbarC = 0.197327053;
const double hbarCsq = hbarC*hbarC;
double KLNModel::Nc = 3.0;
double KLNModel::Nf = 3.0;
double KLNModel::CF = (Nc*Nc-1.0)/(2*Nc);
double KLNModel::Norm=2./CF/hbarCsq;
double KLNModel::Beta0 = (33.0-2.0*Nf)/(12*M_PI);

//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
// This is where the thickness functions are converted into gluon densities
// via the Color-Glass Condensate Model

KLNModel::KLNModel(double srt, int mode, UnintegPartonDist* f)
{
  ecm = srt;
  model = mode;
  wavefunc = f;

  // NN cross sections.
  double sig = hadronxsec::totalXsection(200.0,0);
  double sigel = hadronxsec::elasticXsection(sig,200.0,0,0);
  siginNN200 = sig - sigel;
  sig = hadronxsec::totalXsection(srt,0);
  sigel = hadronxsec::elasticXsection(sig,srt,0,0);
  siginNN = sig - sigel;

  optFixAlpha = 0;   // rc by default
  alphaS = 0.5;      // max value where alpha_s is frozen
  dNdeta = 0;        // default is dN/dy rather than dN/deta
  dEtdy = 0;         // dN/dy is default

  // Parameters for saturation model
  //*******************************************************************
  lambda=0.252;   // evolution speed of Qs(x), for KLN model-ugd

  lambdaQCD=0.2;  // [GeV]
  // HERE
  //lambdaQCD=0.25;  // [GeV]
  lambdaQCD2=lambdaQCD*lambdaQCD;

  mHadron = 0.35;  // mass scale for y->eta conversion
  g2hfac=2.;       // default K factor; may be multiplied by g-->hadr
                   // conversion factor for pt-integr. dN/dy
                   // (not for Et and hadrons from parton fragmentation)

  partonPt = -1.0; // <0 indicates that we integrate over parton pt
  Ptmin=0.05;      // set range of integral over pt (--> dNdy)
  Ptmax=12.0;
}


KLNModel::~KLNModel()
{
}

//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
// Jacobian for dN/dy --> dN/deta transformation; also sets value for y
// corresponding to given eta
double KLNModel::getJacobian(double eta, double pt2)
{
  double msq=mHadron*mHadron;
  double pz=sinh(eta);
  double E=cosh(eta);
  double w=sqrt(msq/pt2+E*E);
  rapidity=0.5*log((w+pz)/(w-pz));
  return E/w;
}


//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
// calculates dNdy or dNdyd2pt for AB collision at rapidity y for proj/targ
// thickness given by npart1, npart2
// this routine is used on the MC sampled nucleon configurations
double KLNModel::getdNdy(double y, double npart1, double npart2, double pt)
{
  Npart1=npart1;
  Npart2=npart2;
  rapidity=y;
  partonPt=pt;  // if >0 then pt is fixed, no integration

  double trEtaY=1.0;
  if(dNdeta<0) {  // global (pt-independent) y->eta Jacobian
    double avg_pt = 0.13 + .32*pow(200./1000.,0.115);
      // this call modifies variable "rapidity" from eta to y !
    trEtaY=getJacobian(rapidity,avg_pt*avg_pt);
  }

  //cout << "->" << ktF_MCintegral() << endl;

  // MC integration over kt(,pt) using bases()
  return trEtaY*g2hfac*Norm*ktF_MCintegral()*9./32.; //Zhi: the 9/32 factor converts the normalization from the pervious version (the line below) to the normalization used in Hirano's CGC code. Details: 9/32=(3/4/pi)^2*(pi^2/2) where the (3/4/pi)^2 is the factor difference of the integrand "func", and the pi^2/2 factor is the difference of the factor in front (g2hfac*Norm).
  // HERE
  // return trEtaY*g2hfac*Norm*ktF_MCintegral();
}


//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
// convolute gluon dN/d2pt with fragm. fct. --> hadron dNd2pt
void KLNModel::FFconv(int iFF, int bins,
		      double ptmin, double dpt, double dNdyd2pt[])
{
  if (iFF<0) {
    cout << "# WARNING: called KLNModel::FFconv() with iFF<0, skipping fragmentation to hadrons!"
	 << endl;
    return;
  }

  int i;
  double* ptbins = NULL;
  ptbins = new double [bins];

  if (ptbins == NULL) {
    cout << "ERROR allocating array in FFconv()" << endl;
    exit(1);
  }
  for (i=0; i<bins; i++)  ptbins[i] = ptmin + i*dpt;

  for (i=0; i<bins; i++) {          // loop over hadron pt
    dNdyd2pt[i] = 0.0;
    double z=1.0, zold=1.0;
    for (int j=i+1; j<bins; j++) {  // integration over parton pt
      z = ptbins[i]/ptbins[j];
      if (z >= 0.05)  // kkp gluon-FF violates mom. sum rule at small z
	dNdyd2pt[i] += (zold-z) * Dgl(iFF,z,ptbins[j]) /z/z * dNdyd2pt[j];
      else break;  // don't need to go to even lower z
      zold = z;
    }
  }
  delete [] ptbins;
}


#include "kkp.cxx"
double KLNModel::Dgl(int iFF, double z, double Q)   // D(z,Q)
{
  int iset = 0;   // LO FF
  int ih = 7;     // charged hadrons
  double dh[11];
  kkpFF(&ih, &iset, &z,	&Q, dh);
  return dh[iFF];  // FF for parton type iFF (0=gluon, see kkp.cxx)
}


//cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
// calls bases MC-integration to obtain dN/dy or dN/deta
double KLNModel::ktF_MCintegral()
{
  bsinit();
  int ndim=3;   // include phi dependence.
  int sample = 5000;

  double x_l[ndim],x_u[ndim];
  int jg[ndim];

  for(int i=0;i<ndim;i++) {
    x_l[i]=0.0;
    x_u[i]=1.0;
    jg[i]=1;
  }

  int nwild=ndim;
  double tune = 1.5;
  int itr1=8;
  int itr2=10;
  double ac1=0.1;
  double ac2=0.1;
  setNoOfSample(sample);
  setTuneValue(tune);
  setIteration1(ac1,itr1);
  setIteration2(ac2,itr2);
  defineVariable(ndim,nwild,x_l,x_u,jg);

  //double x[3] = {0.0, 0.5, 0.5};
  //for (x[0]=0.0; x[0]<1.0; x[0]+=0.1) cout << x[0] << "   " << func(x) << endl;

  // If you want to see the integration step status, set =2
  setPrint(0);
  double aa = bases();
  //double aa = getEstimate();
  //double err = getError();
  return aa;
}


//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
// this function is called by the Bases MC integration routine
// to perform integrals over pt, kt
double KLNModel::func(double* x)
{
  double pt= Ptmin + x[0]*(Ptmax-Ptmin); // integrate over parton pt
  if (partonPt>0.0) pt = partonPt;       // parton pt fixed, no integration
  double ktmax = pt;
  double kt= ktmax*x[1];

  // include phi integral.
  double phi = 2*M_PI*x[2];
  double ktsq1 = 0.25*(pt*pt + kt*kt + 2*kt*pt*cos(phi));
  double ktsq2 = 0.25*(pt*pt + kt*kt - 2*pt*kt*cos(phi));

  transEtaY=1.0;
  double mt=pt;
  double rap_save = rapidity;
  if(dNdeta>0) {
    // this call modifies variable "rapidity" from eta to y !
    transEtaY=getJacobian(rapidity,pt*pt);
    mt=sqrt(pt*pt+mHadron*mHadron);
  }
  double x1=mt/ecm*exp(rapidity);
  double x2=mt/ecm/exp(rapidity);
  rapidity = rap_save;
  if(x1 > 1.0 || x2 > 1.0) return 0.0;

  double qs2a= SaturationScale(x1,Npart1);
  double qs2b= SaturationScale(x2,Npart2);

  double fnc1  = waveFunction(qs2a,x1,ktsq1);
  double fnc2  = waveFunction(qs2b,x2,ktsq2);
  double fnc = fnc1*fnc2;

  //cout << "-->" << qs2a <<"," << x1 << "," << Npart1 << "," << ktsq1 << "," <<fnc1 << endl;

  double scale = max(ktsq1,ktsq2);
  double alp = getAlphaStrong(max(scale,mt*mt));
  //HERE
  //double alp = getAlphaStrong(scale);

  double result = alp*fnc; //alpha_s*phi_A*phi_B
  if (partonPt < 0.0)
    result *= 2.0*M_PI*(Ptmax-Ptmin)*pt; // d^2pT integration
  result /= mt*mt;      // 1/mt^2 factor
  result *= 2.0*M_PI*kt*ktmax;         // d^2kT integration
  result /= 4.0; //Jacobian due to symmetrization of integration
  result *= transEtaY;   // y --> eta Jacobian
  if(dEtdy){
    return result*mt;
  }

  return result;
}

//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
// Qs^2 for _gluon_ probe (adjoint) as a function of thickness
//    thick = # of nucleons / sigma_in(200) * 10  (-> 1/fm^2)
double KLNModel::SaturationScale(double x,double thick) // [GeV^2]
{
  switch (model) {
  case rcBKalbacete :
    // for rcBK: call with Qs_0^2; assumed = 0.4gev^2 for proton at x=0.01
    return thick*siginNN200/10. *0.399;
  case rcBKalbaceteSet2 :
    // for Albacete set 2: Qs_0^2 = 2*0.168gev^2
    return thick*siginNN200/10. *0.336;
  case KLN :
    // KLN: call with Qs(x)^2
    return thick *2./1.53*pow(0.01/x,lambda);
  default :
    cout << "ERROR in KLNModel::SaturationScale(): unknown uGD!\n";
    exit(1);
  }
}


//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
// returns uGD times large-x correction factor
double KLNModel::waveFunction(double qssq, double x, double kt2)
{
  double rapdep=1.0, alp=0.;
  switch (model)
    {
    case KLN:  // KLN model for uGD
      rapdep = pow(1.0-x,4.);
      alp= getAlphaStrong(qssq);
      break;
    case rcBKalbacete: // Albacete's rcBK tabulated ugd, Set 1
    case rcBKalbaceteSet2: // Set 2
      rapdep = pow(1.0-x,4.);
      alp = getAlphaStrong(kt2);
      break;
    default:
      alp = getAlphaStrong(kt2);
    }
  return wavefunc->getFunc(qssq,x,kt2,alp)*rapdep;
}
