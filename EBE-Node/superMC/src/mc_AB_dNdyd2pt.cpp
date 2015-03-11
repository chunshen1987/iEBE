#include "MCnucl.h"
#include "Regge96.h"
#include "OverLap.h"
#include "KLNModel.h"
#include "UnintegPartonDist.h"
#include "KLNfunc.h"


// --- dN/dyd2pt for AB collisions ---

void mc_AB_dNdyd2pt(double ecm, int mode, int lgX, int nevent, int massA, int massB)
{

  // NN cross sections.
  double sig = hadronxsec::totalXsection(ecm,0);
  double sigel = hadronxsec::elasticXsection(sig,ecm,0,0);
  double sigin = sig - sigel;

  cout << "# ecm = " << ecm;
  cout << " mode= " << mode << "   large-x=" << lgX << endl;
  cout << "# A=" << massA << "  B=" << massB << endl;
  cout << "#  sig= " << sig << " sigel= " << sigel << " sigin= " << sigin
       << endl;


  const double maxx=16; // fm  total lattice is [-maxx,maxx]
  const double maxy=16; // fm
  double dx = 0.08; // fm
  double dy = 0.08; // fm

  // set range of impact parameter (note: implement centrality cut via Npart
  // selection; see mc->setCentralityCut() below)
  double bmin=0.;
  double bmax=3.;

  // set no of rapidity bins
  int ny=1;
  double ymin=0., ymax=0.;

  // allocate histo for dN/dyd2pt;
  // ATTN: for rcBK, max pt <~ 20 GeV due to max kt in UGD table phi(x,kt) !
  int iptmax = 200;  // # of pt bins
  double* dNdyd2pt = new double [iptmax];
  for (int ipt=0; ipt<iptmax; ipt++) dNdyd2pt[ipt] =0.0;
  double ptmin = 1.0;   // distribution for pt>ptmin [GeV]
  double dpt = 0.1;     // width of pt bin [GeV]

  int sumNpartA=0;  // N_part(proj) accumulated over events
  int sumNpartB=0;  // N_part(targ) accumulated over events
  int inelEv =0;   // counts inel events
  int sumNcoll=0;  // Ncoll summed over inel events
  
  OverLap* proj = new OverLap(massA,sigin,1); // initialize parameters of
  OverLap* targ = new OverLap(massB,sigin,1); // proj and targ nucleus
  if (massA == 238) {
    proj->setDeformation(1); // U nuclei are deformed
    proj->deformation();
  }
  if (massB == 238) {
    targ->setDeformation(1);
    targ->deformation();
  }


  // initialize proj+targ on transverse grid
  MCnucl* mc = new MCnucl(ecm, maxx, maxy,dx,dy,ny,ymin,ymax);

  // allocate and initialize class for small-x gluons
  KLNModel* kln=0;
  UnintegPartonDist* wf = 0;
  if (mode >= rcBK) wf = new rcBKfunc(mode);
  else wf = new KLNfunc();
  kln = new KLNModel(ecm,mode,wf);

  // modify default K factor (for disc nucleons, rcBKalbacete set 1)
  if (mode == rcBKalbacete) kln->setNormalization(0.75);
  //mc->setNuclProfile(0);      // 0/1: nucleons are Gaussians/discs
  kln->setdNdeta(-1);         // ... global (pt-indep.) Jacobian
  mc->setKLN(kln);

  // set range of accepted Npart (centrality cut)
  mc->setCentralityCut(10,420);

  // generates look-up table of dN/dyd2pt as fct of proj/targ thicknesses
  mc->makeTable(ptmin,dpt,iptmax);

  // start event loop
  cout << "#\n# starting event loop\n";
  for(int iev=0;iev<nevent; iev++) {
    double b = sqrt((bmax*bmax - bmin*bmin)*drand48()+bmin*bmin);
    mc->generateNucleus(b,proj,targ);
    int binary = mc->getBinaryCollision();
    if(binary==0 || mc->CentralityCut()==0) {
      mc->deleteNucleus();
      continue;  // no interaction occurs.
    }
    inelEv++;
    sumNpartA += mc->getNpart1();
    sumNpartB += mc->getNpart2();
    sumNcoll += binary;

    mc->getTA2(); // init density of "valence" charges per unit transv. area
    for (int ipt=0; ipt<iptmax; ipt++) {
      mc->setDensity(0,ipt);   // rapidity and pt bin
      // dN/dyd2pt via \int d^2rt, summed over events
      dNdyd2pt[ipt] += mc->getdNdy();
    }
    mc->deleteNucleus();
  } // event loop

  int iFF=-1;    // fragm. fct: -1=off, 0=gluon
  double Etparton=0.;   // dEt/dy before fragmentation
  double Ethadron=0.;   // dEt/dy after fragmentation
  for (int ipt=0; ipt<iptmax; ipt++) {
    Etparton += 2.*M_PI*dpt*pow(ptmin+dpt*ipt,2.)*dNdyd2pt[ipt];
  }

  // convolute gluon spectrum with fragm. fct.
  if (iFF >= 0)
    kln->FFconv(iFF,iptmax,ptmin,dpt,dNdyd2pt);
  for (int ipt=0; ipt<iptmax; ipt++)
    Ethadron += 2.*M_PI*dpt*pow(ptmin+dpt*ipt,2.)*dNdyd2pt[ipt];

  cout << "# " << inelEv << " inel. events\n#\n";
  cout << "# <Npart_A>   <Npart_B>   <Ncoll>: "<< 
    (double)sumNpartA/inelEv << " " <<
    (double)sumNpartB/inelEv << " " <<
    (double)sumNcoll/inelEv  << endl;

  if (inelEv>0) {
    cout << "# FF for parton type " << iFF << endl;
    cout << "# Et(g) = " << Etparton/inelEv
	 << "    Et(h)/Et(g) = " << Ethadron/Etparton << endl;
    cout << "# pt\tdNch/d2pt" << endl;
    for (int ipt=0; ipt<iptmax; ipt++)
      cout << setprecision(3) << setw(5) << ipt*dpt+ptmin
	   << setprecision(4) << setw(11) << dNdyd2pt[ipt]/inelEv
	   << endl;
  }
  delete [] dNdyd2pt;
  delete proj; delete targ;
  delete kln;  delete mc;  delete wf;
}
