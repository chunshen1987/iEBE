#include "ParamDefs.h"
#include "MCnucl.h"
#include "Regge96.h"
#include "OverLap.h"
#include "KLNModel.h"
#include "UnintegPartonDist.h"
#include "KLNfunc.h"



// --- dN/dyd2pt for pA collision from kt-factorization ---

void ktFpA(double ecm, int mode, int lgX, int nevent, int massA)
{

  // NN cross sections.
  double sig = hadronxsec::totalXsection(ecm,0);
  double sigel = hadronxsec::elasticXsection(sig,ecm,0,0);
  double sigin = sig - sigel;

  cout << "#pA,  sqrt_s= " << ecm << endl;
  cout << "# A=" << massA << endl;
  cout << "#sig_in= " << sigin << endl;
  cout << "# mode= " << mode << "   large-x=" << lgX << endl;

  const double maxx=8.; // fm  total lattice is [-maxx,maxx]
  const double maxy=8.; // fm
  double dx = 0.08; // fm
  double dy = 0.08; // fm

  // set range of impact parameter
  double bmin=0.;
  double bmax=8.;
  cout << "# b_max = " << bmax << endl;

  int sumNcoll=0;  // Ncoll summed over inel events
  int sumNpartA=0;  // N_part(proj) accumulated over events
  int sumNpartB=0;  // N_part(targ) accumulated over events
  int inelEv =0;      // counts inel events

  OverLap* proj = new OverLap(massP,sigin); // initialize parameters of
  OverLap* targ = new OverLap(massA,sigin); // proj and targ nucleus

  // initialize proj+targ on transverse grid
  int iy=1;             // # of rapidity bins
  double ymin = 3.2;    // min rapidity
  double ymax = 3.2;    // max rapidity
  MCnucl* mc = new MCnucl(ecm, maxx, maxy,dx,dy,iy,ymin,ymax);
  mc->setNuclProfile(0);      // 0/1: nucleons are Gaussians/discs

  KLNModel* kln=0;   // class for small-x partons
  UnintegPartonDist* wf = 0;  // class for uGD
  if (mode >= rcBK) wf = new rcBKfunc(mode);
  else wf = new KLNfunc();
  kln = new KLNModel(ecm,mode,wf);

  //kln->setNormalization(0.75);  // for disc shape nucleons, rcBKalbacete uGD
  //kln->setdNdeta(1);   // flag for eta instead of y; local (pt-dep.) Jacobian
  kln->setdNdeta(-1);         // ... global (pt-indep.) Jacobian
  mc->setKLN(kln);

  // allocate histo for dN/dyd2pt;
  // ATTN: for rcBK, max pt <~ 20 GeV due to max kt in UGD table phi(x,kt) !
  int iptmax = 200;  // # of pt bins
  double* dNdyd2pt = new double [iptmax];
  for (int ipt=0; ipt<iptmax; ipt++) dNdyd2pt[ipt] =0.0;
  double ptmin = 1.0;   // distribution for pt>ptmin [GeV]
  double dpt = 0.1;     // width of pt bin [GeV]

  kln->setPtRange(ptmin, ptmin+iptmax*dpt);
  Large_x* val = 0;
  if (lgX && mode>=rcBK) {    // large-xF 'valence' partons (rcBK mode only)
    val = new Large_x(kln,wf);  mc->setLgX(val);
  }

  // generates look-up table of dN/dyd2ptd2rt as fct of proj/targ thicknesses
  mc->makeTable(ptmin,dpt,iptmax);

  // start event loop
  cout << "#\n# starting event loop\n";
  for (int iev=0; iev<nevent; iev++) {
    double b = sqrt((bmax*bmax - bmin*bmin)*drand48()+bmin*bmin);
    mc->generateNucleus(b,proj,targ);
    int binary = mc->getBinaryCollision();
    if (binary != 0) {    // if interaction occurs
      inelEv++;
      //cout << "# event " << inelEv << "  " << "  b=" << b << endl;
      sumNpartA += mc->getNpart1();
      sumNpartB += mc->getNpart2();
      sumNcoll += binary;
      mc->getTA2(); // init density of "valence" charges per unit transv. area
      for (int ipt=0; ipt<iptmax; ipt++) {
	mc->setDensity(0,ipt);  // rapidity bin and pt bin
	// dN/dyd2pt via \int d^2rt, summed over events
	dNdyd2pt[ipt] += mc->getdNdy();
      }
    }
    mc->deleteNucleus();
  } // event loop


  int iFF=0;    // fragm. fct: -1=off, 0=gluon
  if (lgX) iFF = -100;  // ATTN: -lgX flag already activated fragmentation !
  double Etparton=0.;   // dEt/dy before fragmentation
  double Ethadron=0.;   // dEt/dy after fragmentation
  for (int ipt=0; ipt<iptmax; ipt++) {
    Etparton += 2.*M_PI*dpt*pow(ptmin+dpt*ipt,2.)*dNdyd2pt[ipt];
  }

  // convolute small-x gluon spectrum with fragm. fct.
  if (iFF >= 0) {
    kln->FFconv(iFF,iptmax,ptmin,dpt,dNdyd2pt);
    for (int ipt=0; ipt<iptmax; ipt++)
      Ethadron += 2.*M_PI*dpt*pow(ptmin+dpt*ipt,2.)*dNdyd2pt[ipt];
  }

  cout << "# " << inelEv << " inel. events\n#\n";
  if (inelEv>0) {
    cout << "# FF for parton type " << iFF << endl;
    if (iFF >= 0) cout << "# Et(g) = " << Etparton/inelEv
	 << "    Et(h)/Et(g) = " << Ethadron/Etparton << endl;

    // output: <Npart_A> <Npart_B> <Ncoll> over all events
    cout << "# <Npart_A>  <Npart_B>  <Ncoll>" << endl;
    cout << "# " << (double)sumNpartA/inelEv << " " <<
      (double)sumNpartB/inelEv << " " <<
      (double)sumNcoll/inelEv << endl;

    cout << "#\n# pt\tdNch/d2pt" << endl;
    for (int ipt=0; ipt<iptmax; ipt++)
      cout << setprecision(3) << setw(5) << ipt*dpt+ptmin
	   << setprecision(4) << setw(11) << dNdyd2pt[ipt]/inelEv
	   << endl;
  }
  delete [] dNdyd2pt;
  delete proj; delete targ;
  delete kln;  delete mc;  delete wf;
}
