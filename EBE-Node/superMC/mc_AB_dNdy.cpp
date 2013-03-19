#include "MCnucl.h"
#include "Regge96.h"
#include "OverLap.h"
#include "KLNModel.h"
#include "UnintegPartonDist.h"
#include "KLNfunc.h"


// --- dN/dy (or dEt/dy, eccentricity) for AB collisions ---

void mc_AB_dNdy(double ecm, int mode, int lgX, int Et, int nevent, int massA, int massB)
{

  // NN cross sections.
  double sig = hadronxsec::totalXsection(ecm,0);
  double sigel = hadronxsec::elasticXsection(sig,ecm,0,0);
  double sigin = sig - sigel;

  cout << "# ecm = " << ecm;
  cout << " mode= " << mode << "   large-x=" << lgX << endl;
  cout << "# A=" << massA << "  B=" << massB << endl;
  cout << "# tr. energy flag Et=" << Et << endl;
  cout << "#  sig= " << sig << " sigel= " << sigel << " sigin= " << sigin
       << endl;


  const double maxx=16; // fm  total lattice is [-maxx,maxx]
  const double maxy=16; // fm
  double dx = 0.08; // fm
  double dy = 0.08; // fm

  // set range of impact parameter (note: implement centrality cut via Npart
  // selection; see mc->setCentralityCut() below)
  double bmin=0.;
  double bmax=14.;

  // eccentricity of produced gluons from dN/dyd2rt
  double eccPart, dummy; 

  // set no of rapidity bins
  int ny=1;
  double ymin=0., ymax=0.;

  double *dndy = new double [ny];
  double *sum_dndy = new double [ny];
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
  kln->setdEtdy(Et);         // compute Et or N
  // g-->hadr multiplication factor for pt-integr. dN/dy, LHC energy
  if (Et==0) kln->setg2hfac(5.0);
  //kln->setdNdeta(1);    // flag for eta instead of y; local (pt-dep.) Jacobian
  kln->setdNdeta(-1);         // ... global (pt-indep.) Jacobian
  mc->setKLN(kln);

  Large_x* val = 0;
  if (lgX && mode>=rcBK) {    // large-xF partons (rcBK mode only)
    val = new Large_x(kln,wf);  mc->setLgX(val);
  }

  // set range of accepted Npart (centrality cut)
  mc->setCentralityCut(10,420);

  // generates look-up table of dN/dy as fct of proj/targ thicknesses
  mc->makeTable();

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
    for(int iy=0;iy<ny;iy++) { // loop over rapidity bins
      // initializes dN/dyd2rt on 2d transv. grid (for this event)
      // can be accessed via mc->getRho(iy, x, y)
      mc->setDensity(iy);
      dndy[iy] = mc->getdNdy(); // dN/dy via \int d^2rt
      sum_dndy[iy] += dndy[iy];
      // eccentricity of produced gluons at y=0
      mc->getEccentricityGrid(dummy, eccPart, dummy, dummy, iy);

      cout 
	<< setprecision(3) << setw(5) << b
	<< setprecision(3) << setw(7) << mc->getNpart1()+mc->getNpart2()
	<< setprecision(3) << setw(6) << ymin+(ymax-ymin)/ny*iy
	<< setprecision(4) << setw(10)<< dndy[iy]
	<< setprecision(3) << setw(9)<< eccPart
	<< endl;
    } // end loop over rapidity bins
    mc->deleteNucleus();
  } // event loop

  // output: <Npart_A> <Npart_B> <Ncoll> <dN/dy> over all events
  for(int iy=0;iy<ny;iy++)
    cout << "# " << (double)sumNpartA/inelEv << " " <<
      (double)sumNpartB/inelEv << " " <<
      (double)sumNcoll/inelEv << " " << ymin+(ymax-ymin)/ny*iy
	 << " " << sum_dndy[iy]/inelEv << endl;

  if(wf) delete wf;
  if(kln) delete kln;
  if (val) delete val;
  delete proj; delete targ;
  delete mc;
}
