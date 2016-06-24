#include "MCnucl.h"
#include "Regge96.h"
#include <iomanip>
#include <ctime>
#include <sys/time.h>

#include "arsenal.h"
#include "GaussianNucleonsCal.h"
#include "ParameterReader.h"
#include "NBD.h"

#define HBARC 0.197327053

using namespace std;

//*********************************************************************
// functions to initialize the 2d transv. coordinate grid,
// place two colliding nuclei on the grid, determine their thickness
// generate lookup table for dN/dy as fct of T_A(rt), T_B(rt)

MCnucl::MCnucl(ParameterReader* paraRdr_in) {
    paraRdr = paraRdr_in;

    // default Alpha is Glauber model: sd = (1-Alpha)WN + (Alpha)BC
    Alpha = paraRdr->getVal("alpha");

    // tmax-1 is max # of locally overlapping nucleons,
    // upper limit in dNdy table
    tmax = paraRdr->getVal("tmax");

    // fix grid properties
    Xmax = paraRdr->getVal("maxx");
    Ymax = paraRdr->getVal("maxy");
    Xmin = -Xmax;
    Ymin = -Ymax;
    dx = paraRdr->getVal("dx");
    dy = paraRdr->getVal("dy");
    Maxx=(int)((Xmax-Xmin)/dx+0.1)+1;
    Maxy=(int)((Ymax-Ymin)/dy+0.1)+1;

    // default npart cutoff
    NpartMax = paraRdr->getVal("Npmax");
    NpartMin = paraRdr->getVal("Npmin");

    // Unintegrated PT optns
    PTinte = paraRdr->getVal("PT_flag");
    PTmax  = paraRdr->getVal("PT_Max");
    PTmin  = paraRdr->getVal("PT_Min");
    dpt = paraRdr->getVal("d_PT");
    MaxPT=(int)((PTmax-PTmin)/dpt+0.1)+1;
    if (PTinte < 0)
        PT_order = paraRdr->getVal("PT_order");   
    else
        PT_order = 1;  // does not apply when there is no PT integration
  
    //.... NN cross sections in mb
    double ecm = paraRdr->getVal("ecm");
    double sig = hadronxsec::totalXsection(200.0,0);
    double sigel = hadronxsec::elasticXsection(sig,200.0,0,0);
    siginNN200 = sig - sigel;
    sig = hadronxsec::totalXsection(ecm,0);
    sigel = hadronxsec::elasticXsection(sig,ecm,0,0);
    siginNN = sig - sigel;

    // Include additional CC fluctution?
    CCFluctuationModel = paraRdr->getVal("cc_fluctuation_model");
    CCFluctuationK = paraRdr->getVal("cc_fluctuation_k");
    if (CCFluctuationModel) nbd = new NBD;
    if (CCFluctuationModel > 5) {
        gsl_rng_env_setup();
        gslRngType = gsl_rng_default;
        gslRng = gsl_rng_alloc(gslRngType);
        timeval a;
        gettimeofday(&a, 0);
        // randomSeed use CPU clock
        int randomSeed=a.tv_usec;
        //initialize random generator
        gsl_rng_set (gslRng, (unsigned long int) randomSeed);
        ccFluctuationGammaTheta = paraRdr->getVal("cc_fluctuation_Gamma_theta");
    }

    which_mc_model = paraRdr->getVal("which_mc_model");
    sub_model = paraRdr->getVal("sub_model");
    shape_of_nucleons = paraRdr->getVal("shape_of_nucleons");

    gaussCal = NULL;
    entropy_gaussian_width = 0.0;
    paraRdr->setVal("siginNN", siginNN);
    // for Gaussian-shaped nucleons calculations
    gaussCal = new GaussianNucleonsCal(paraRdr);
    entropy_gaussian_width = gaussCal->entropy_gaussian_width;
    entropy_gaussian_width_sq = entropy_gaussian_width*entropy_gaussian_width;

    // adding quark substructure Fluctuations (from Kevin Welsh)
    shape_of_entropy = paraRdr->getVal("shape_of_entropy");
    // For separation of entropy to collision detection (Kevin)
    if (shape_of_entropy == 3) {
        quark_width = paraRdr->getVal("quark_width");
        double nucleon_width = gaussCal->width;
        quark_dist_width = sqrt(nucleon_width*nucleon_width - quark_width*quark_width);
        gaussDist = new GaussianDistribution(0, quark_dist_width);
    }

    dndyTable=0;    // lookup table pointers not valid yet
    dndydptTable=0;
    overSample=1;  // default: no oversampling
    binRapidity = paraRdr->getVal("ny");
    rapMin = paraRdr->getVal("rapMin");
    rapMax = paraRdr->getVal("rapMax");
    rapidity = 0.0;

    dsq = 0.1*siginNN/M_PI/overSample;

    kln=0;  // pointer to class for small-x gluons
    val=0;  // pointer to class for large-x (x>x0) partons

    TA1 = new double* [Maxx];    // 2d grid for proj/targ. thickness functions
    TA2 = new double* [Maxx];
    sd_TA1 = new double* [Maxx];
    sd_TA2 = new double* [Maxx];
    rho_binary = new double* [Maxx];  // 2d grid for the binary density
    // 2d grid for spectator density from nucleus 1
    spectator_1 = new double* [Maxx];
    // 2d grid for spectator density from nucleus 2
    spectator_2 = new double* [Maxx];
    for (int ix = 0; ix < Maxx; ix++) {
        TA1[ix] = new double[Maxy];
        TA2[ix] = new double[Maxy];
        sd_TA1[ix] = new double[Maxy];
        sd_TA2[ix] = new double[Maxy];
        rho_binary[ix] = new double[Maxy];
        spectator_1[ix] = new double[Maxy];
        spectator_2[ix] = new double[Maxy];
        for (int iy = 0; iy < Maxy; iy++) {
            TA1[ix][iy] = 0;
            TA2[ix][iy] = 0;
            sd_TA1[ix][iy] = 0.;
            sd_TA2[ix][iy] = 0.;
            rho_binary[ix][iy] = 0;
            spectator_1[ix][iy] = 0.;
            spectator_2[ix][iy] = 0.;
        }
    }

    rho = new GlueDensity(Xmax, Ymax, PTmin, PTmax,
                          dx, dy, dpt, binRapidity, rapMin,rapMax);

    Xcm=0.0, Ycm=0.0, angPart=0.0;

    // flag to include nucleon-nucleon correlation 
    // when sampled the nucleon positions
    flag_NN_correlation = paraRdr->getVal("include_NN_correlation");
}


MCnucl::~MCnucl() {
    for (int ix = 0; ix < Maxx; ix++) {
        delete [] TA1[ix];
        delete [] TA2[ix];
        delete [] sd_TA1[ix];
        delete [] sd_TA2[ix];
        delete [] rho_binary[ix];
        delete [] spectator_1[ix];
        delete [] spectator_2[ix];
    }
    delete [] TA1;
    delete [] TA2;
    delete [] sd_TA1;
    delete [] sd_TA2;
    delete [] rho_binary;
    delete [] spectator_1;
    delete [] spectator_2;

    delete rho;

    if (dndyTable) {
        for (int iy = 0; iy < binRapidity; iy++) {
            for(int j = 0; j < tmax; j++)
                delete [] dndyTable[iy][j];
            delete [] dndyTable[iy];
        }
        delete [] dndyTable;
    }

    if (dndydptTable) {
        for (int iy = 0; iy < binRapidity; iy++) {
            for (int j = 0; j < tmaxPt; j++) {
                for (int i = 0; i < tmaxPt; i++)
                    delete [] dndydptTable[iy][j][i];
                delete [] dndydptTable[iy][j];
            }
            delete [] dndydptTable[iy];
        }
        delete [] dndydptTable;
    }

    if (CCFluctuationModel > 5) gsl_rng_free(gslRng);

    if (gaussCal) delete gaussCal;
}



/* place two nuclei/nucleons on the transverse lattice (separated by
   impact parameter b)
*/
void MCnucl::generateNucleus(double b, OverLap* proj, OverLap* targ)
{
  double rmin=0.9*0.9; // minimal nucleon separation (squared; in fm^2).
  double cx, phi;
  for(int ie=0;ie<overSample;ie++) {

    // nucleus 1
    double xcm=0.0, ycm=0.0, zcm=0.0;
    int nn=proj->getAtomic();
    cx = 1.0-2.0*drand48();
    phi = 2*M_PI*drand48();
    lastCx1=cx;
    lastPh1=phi;
    proj ->setRotation(cx, phi);
    if (nn == 1)
    {
      nucl1.push_back(new Particle(b/2.0, 0.0, 0.0)); // shift along x-axis
    }
    else if (nn==2) // in the case of a deuteron (added by Brian Baker)
    {
      double x1, y1, z1, x2, y2, z2;
         
      proj-> GetDeuteronPosition(x1,y1,z1,x2,y2,z2);
      nucl1.push_back(new Particle(x1+(b/2.0),y1,z1)); // shift along x-axis
      nucl1.push_back(new Particle(x2+(b/2.0),y2,z2));
    }
    else if(nn==3) // in the case of 3He
    {
      double x1, x2, x3, y1, y2, y3, z1, z2, z3;
      proj->GetTritonPosition(x1,y1,z1,x2,y2,z2,x3,y3,z3); // Triton data points from Carlson have center of mass (0,0,0).
      nucl1.push_back(new Particle(x1+(b/2.0),y1,z1));
      nucl1.push_back(new Particle(x2+(b/2.0),y2,z2));
      nucl1.push_back(new Particle(x3+(b/2.0),y3,z3));
    }
    else 
    {
       if(flag_NN_correlation == 0 || nn != 197 || nn != 208)
       {
          // for large nuclei
          for(int ia=0;ia<nn;ia++) {
            double x,y,z;
            int icon=0;
            do {
              proj->getDeformRandomWS(x,y,z);
              icon=0;
              for(int i=ie*nn; i<(int)nucl1.size();i++) {
                double x1=nucl1[i]->getX();
                double y1=nucl1[i]->getY();
                double z1=nucl1[i]->getZ();
                double r2 = (x-x1)*(x-x1) + (y-y1)*(y-y1) + (z-z1)*(z-z1);
                if(r2 < rmin) {
                  icon=1;
                  break;
                }
              }
            } while(icon==1);

            xcm +=x;
            ycm +=y;
            zcm +=z;
            nucl1.push_back(new Particle(x,y,z));
          }

          for(int ia=0;ia<nn;ia++) { // shift center of nucleus
            double x = nucl1[ia]->getX() - xcm/nn + b/2.0;
            double y = nucl1[ia]->getY() - ycm/nn;
            double z = nucl1[ia]->getZ() - zcm/nn;
            nucl1[ia]->setX(x);
            nucl1[ia]->setY(y);
            nucl1[ia]->setZ(z);
          }
       }
       else  // load nucleon positions from file
       {
          double **nucleon_positions = new double* [nn];
          for(int inucleon = 0; inucleon < nn; inucleon++)
          {
             nucleon_positions[inucleon] = new double [3];
          }
          proj->get_nucleon_position_with_NN_correlation(nucleon_positions);

          double xcm = 0.0, ycm = 0.0, zcm = 0.0;
          for(int ia = 0; ia < nn; ia++)
          {
             double x_local = nucleon_positions[ia][0];
             double y_local = nucleon_positions[ia][1];
             double z_local = nucleon_positions[ia][2];
             xcm += x_local;
             ycm += y_local;
             zcm += z_local;
             nucl1.push_back(new Particle(x_local, y_local, z_local));
          }
          
          for(int ia = 0; ia < nn; ia++)
          {
             double x_local = nucl1[ia]->getX() - xcm/nn + b/2.0;
             double y_local = nucl1[ia]->getY() - ycm/nn;
             double z_local = nucl1[ia]->getZ() - zcm/nn;
             nucl1[ia]->setX(x_local);
             nucl1[ia]->setY(y_local);
             nucl1[ia]->setZ(z_local);
          }

          // clean up
          for(int inucleon = 0; inucleon < nn; inucleon++)
          {
             delete [] nucleon_positions[inucleon];
          }
          delete [] nucleon_positions;
       }
    }


    // nucleus 2
    xcm=0.0; ycm=0.0; zcm=0.0;
    nn=targ->getAtomic();

    cx = 1.0-2.0*drand48();
    phi = 2*M_PI*drand48();
    targ->setRotation(cx, phi);
    lastCx2=cx;
    lastPh2=phi;
    if (nn == 1)
    {
      nucl2.push_back(new Particle(-b/2.0, 0.0, 0.0)); // shift along x-axis
    }
    else if (nn==2) // in the case of a deuteron (added by Brian Baker)
    {
      double x1, y1, z1, x2, y2, z2;
         
      targ-> GetDeuteronPosition(x1,y1,z1,x2,y2,z2);
      nucl2.push_back(new Particle(x1-(b/2.0),y1,z1)); // shift along x-axis
      nucl2.push_back(new Particle(x2-(b/2.0),y2,z2));
    }
    else if(nn==3) // in the case of 3He
    {
      double x1, x2, x3, y1, y2, y3, z1, z2, z3;
      targ->GetTritonPosition(x1,y1,z1,x2,y2,z2,x3,y3,z3); // Triton data points from Carlson have center of mass (0,0,0).
      nucl2.push_back(new Particle(x1-(b/2.0),y1,z1));
      nucl2.push_back(new Particle(x2-(b/2.0),y2,z2));
      nucl2.push_back(new Particle(x3-(b/2.0),y3,z3));
    }
    else 
    {
       if(flag_NN_correlation == 0 || nn != 197 || nn != 208)
       {
          for(int ia=0;ia<nn;ia++) {
            double x,y,z;
            int icon=0;
            do {
              targ->getDeformRandomWS(x,y,z);
              icon=0;
              for(int i=ie*nn; i<(int)nucl2.size();i++) {
                double x1=nucl2[i]->getX();
                double y1=nucl2[i]->getY();
                double z1=nucl2[i]->getZ();
                double r2 = (x-x1)*(x-x1) + (y-y1)*(y-y1) + (z-z1)*(z-z1);
                if(r2 < rmin) {
                  icon=1;
                  break;
                }
              }
            } while(icon==1);

            xcm +=x;
            ycm +=y;
            zcm +=z;
            nucl2.push_back(new Particle(x,y,z));
          }

          for(int ia=0;ia<nn;ia++) {
            double x = nucl2[ia]->getX() - xcm/nn - b/2.0;
            double y = nucl2[ia]->getY() - ycm/nn;
            double z = nucl2[ia]->getZ() - zcm/nn;
            nucl2[ia]->setX(x);
            nucl2[ia]->setY(y);
            nucl2[ia]->setZ(z);
          }
       }
       else  // load nucleon positions from file
       {
          double **nucleon_positions = new double* [nn];
          for(int inucleon = 0; inucleon < nn; inucleon++)
          {
             nucleon_positions[inucleon] = new double [3];
          }
          targ->get_nucleon_position_with_NN_correlation(nucleon_positions);

          double xcm = 0.0, ycm = 0.0, zcm = 0.0;
          for(int ia = 0; ia < nn; ia++)
          {
             double x_local = nucleon_positions[ia][0];
             double y_local = nucleon_positions[ia][1];
             double z_local = nucleon_positions[ia][2];
             xcm += x_local;
             ycm += y_local;
             zcm += z_local;
             nucl2.push_back(new Particle(x_local, y_local, z_local));
          }

          for(int ia = 0; ia < nn; ia++)
          {
             double x_local = nucl2[ia]->getX() - xcm/nn - b/2.0;
             double y_local = nucl2[ia]->getY() - ycm/nn;
             double z_local = nucl2[ia]->getZ() - zcm/nn;
             nucl2[ia]->setX(x_local);
             nucl2[ia]->setY(y_local);
             nucl2[ia]->setZ(z_local);
          }

          // clean up
          for(int inucleon = 0; inucleon < nn; inucleon++)
          {
             delete [] nucleon_positions[inucleon];
          }
          delete [] nucleon_positions;
       }
    }
  }

}



// --- find participants from proj/target and the number of binary coll. ---
int MCnucl::getBinaryCollision()
{
  int* mapping_table1 = new int[(int)nucl1.size()]; // it stores the index of the nucleon in the participant array if it is wounded; otherwise it's set to -1 (next line)
  for (int i=0; i<(int)nucl1.size(); i++) mapping_table1[i]=-1; // -1 means it's not wounded
  int* mapping_table2 = new int[(int)nucl2.size()];
  for (int i=0; i<(int)nucl2.size(); i++) mapping_table2[i]=-1; // -1 means it's not wounded
  // decide collision pairs
  Ncoll=0;
  for(int i=0;i<(int)nucl1.size();i++) { // loop over proj. nucleons
    double x1 = nucl1[i]->getX();
    double y1 = nucl1[i]->getY();
    for(int j=0;j<(int)nucl2.size();j++) { // loop over targ. nucleons
      double x2 = nucl2[j]->getX();
      double y2 = nucl2[j]->getY();
      double dc = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
      if(hit(sqrt(dc))) {
        Ncoll++;
        // Take care of wounded nucleons registration:
        nucl1[i]->setNumberOfCollision();
        // push targ. nucleon j to participant stack (only once though)
        if(nucl1[i]->getNumberOfCollision()==1)
        {
          participant.push_back(new Participant(nucl1[i],1));
          if(CCFluctuationModel > 5)
          {
            if(shape_of_entropy == 3)
            {
              double f1 = sampleFluctuationFactorforParticipant();
              double f2 = sampleFluctuationFactorforParticipant();
              double f3 = sampleFluctuationFactorforParticipant();
              participant.back()->getParticle()->setfluctfactorQuarks(f1,f2,f3);
            }
            else
              participant.back()->setfluctfactor(sampleFluctuationFactorforParticipant());
          }
          mapping_table1[i] = participant.size()-1;
        }
        nucl2[j]->setNumberOfCollision();
        // push targ. nucleon j to participant stack (only once though)
        if(nucl2[j]->getNumberOfCollision()==1)
        {
          participant.push_back(new Participant(nucl2[j],2));
          if(CCFluctuationModel > 5)
          {
            if(shape_of_entropy == 3)
            {
              double f1 = sampleFluctuationFactorforParticipant();
              double f2 = sampleFluctuationFactorforParticipant();
              double f3 = sampleFluctuationFactorforParticipant();
              participant.back()->getParticle()->setfluctfactorQuarks(f1,f2,f3);
            }
            else
              participant.back()->setfluctfactor(sampleFluctuationFactorforParticipant());
          }
          mapping_table2[j] = participant.size()-1;
        }
        // Take care of binary collision registration:
        binaryCollision.push_back(new CollisionPair((x1+x2)/2,(y1+y2)/2));
        if(CCFluctuationModel > 5)
           binaryCollision.back()->setfluctfactor(sampleFluctuationFactorforBinaryCollision());
        if (which_mc_model==5 && sub_model==2) // need to know which binary collision happened to which participants
        {
          int current_binaryCollision_index = binaryCollision.size()-1;
          // record collisions to the current participants
          participant[mapping_table1[i]]->who_hit_me.push_back(current_binaryCollision_index);
          participant[mapping_table2[j]]->who_hit_me.push_back(current_binaryCollision_index);
        }
      }
    }
  }

  // Additional treatment for Uli-Glb model; restore collision info as weight into binaryCollision array
  if (which_mc_model==5 && sub_model==2) // Uli-Glb model
  {
    for (int i=0; i<participant.size(); i++) // loop over wounded nucleons
    {
      int number_of_collsions_to_this_nucleon = participant[i]->who_hit_me.size();
      double weight_to_add = 1./number_of_collsions_to_this_nucleon;
      for (int j=0; j<number_of_collsions_to_this_nucleon; j++)
      {
        binaryCollision[participant[i]->who_hit_me[j]]->additional_weight += weight_to_add;
      }
    } // <-> for (int i=0; i<participant.size(); i++)
  } // <-> (which_mc_model==5 && sub_model==2)

  int npart = participant.size();
  Npart1=0;
  Npart2=0;
  for(int i=0;i<npart;i++) {
    if(participant[i]->isNucl() == 1) Npart1++;
    if(participant[i]->isNucl() == 2) Npart2++;
  }
  if(npart != Npart1+Npart2) {
    cout << " something is wrong with # of participants " << endl;
    exit(1);
  }
  //cout << "Npart1=" << Npart1 << ", " << "Npart2=" << Npart2 << endl;
  Npart1 /= overSample;
  Npart2 /= overSample;

  delete[] mapping_table1;
  delete[] mapping_table2;

  return binaryCollision.size();
}

/* Original
int MCnucl::getBinaryCollision()
{
  Ncoll=0;
  for(int i=0;i<(int)nucl1.size();i++) { // loop over proj. nucleons
    double x1 = nucl1[i]->getX();
    double y1 = nucl1[i]->getY();
    int part=0;  // no of collision partners for proj. nucleon i
    for(int j=0;j<(int)nucl2.size();j++) { // loop over targ. nucleons
      double x2 = nucl2[j]->getX();
      double y2 = nucl2[j]->getY();
      double dc = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
      if(hit(sqrt(dc))) {
        part++;
        Ncoll++;
        binaryCollision.push_back(new CollisionPair((x1+x2)/2,(y1+y2)/2));
        nucl2[j]->setNumberOfCollision();
        // push targ. nucleon j to participant stack (only once though)
        if(nucl2[j]->getNumberOfCollision()==1)
          participant.push_back(new Participant(nucl2[j],2));
      }
    }
    if(part) { // push proj. nucleon i to participant stack, if needed
      nucl1[i]->setNumberOfCollision(part);
      participant.push_back(new Participant(nucl1[i],1));
    }
  }

  int npart = participant.size();
  Npart1=0;
  Npart2=0;
  for(int i=0;i<npart;i++) {
    if(participant[i]->isNucl() == 1) Npart1++;
    if(participant[i]->isNucl() == 2) Npart2++;
  }
  if(npart != Npart1+Npart2) {
    cout << " something is wrong with # of participants " << endl;
    exit(1);
  }
  //cout << "Npart1=" << Npart1 << ", " << "Npart2=" << Npart2 << endl;
  Npart1 /= overSample;
  Npart2 /= overSample;

  return binaryCollision.size();
}*/



// old stuff
//  (YN): determine whether nucleons separated by distance dr2 interact or not
// interaction probability at impact param. b is 1-exp(-sigeff(s)*Tpp(b))
int MCnucl::hit(double b)
{
  if (shape_of_nucleons==1) // disc
    return (b*b<=dsq) ? 1 : 0;  // |r1-r2|^2 < sigma_in(s) / pi
  else if (shape_of_nucleons>=2 && shape_of_nucleons<=9) // Gaussian nucleon shape
    return gaussCal->testCollision(b);
  else
  {
    cout << "MCnucl::hit: you shouldn't come to this line. Check your parameters" << endl;
    exit(-1);
  }
}



// checks whether Npart1+Npart2 is in the desired range
//   ATTN: Npart1, Npart2 need to be initialized through getBinaryCollision() !
int MCnucl::CentralityCut()
{
  int Nptot = Npart1 + Npart2;
  if (Nptot<=NpartMax && Nptot>=NpartMin) return 1;
  return 0;
}


// --- determine thickness of proj+targ nuclei over 2d transv. grid,
//     for given MC event ---
void MCnucl::getTA2() {
    int npart = participant.size();
    double areai = 10.0/siginNN;
    int imax = 0;

    double nucleon_width;
    if (shape_of_nucleons >= 2 && shape_of_nucleons <= 9)
        nucleon_width = gaussCal->width;
    double dc_sq_max_gaussian = 25.*nucleon_width*nucleon_width;
    double d_max;
    if (shape_of_nucleons == 1)
        d_max = 2.*sqrt(dsq);
    if (shape_of_nucleons >= 2 && shape_of_nucleons <=9)
        d_max = 5.*nucleon_width;
    double *part_x, *part_y;
    part_x = new double [npart];
    part_y = new double [npart];
    for (int i = 0; i < npart; i++) {
        part_x[i] = participant[i]->getX();
        part_y[i] = participant[i]->getY();
    }

    for (int ix = 0; ix < Maxx; ix++) {
        for (int iy = 0; iy < Maxy; iy++) {
            TA1[ix][iy] = 0.0;
            TA2[ix][iy] = 0.0;
        }
    }

    // loop over nucleons which overlap the grid point (xg,yg)
    for (unsigned int ipart=0; ipart<npart; ipart++) {
        double x = part_x[ipart];
        double y = part_y[ipart];
        int x_idx_left = (int)((x - d_max - Xmin)/dx);
        int x_idx_right = (int)((x + d_max - Xmin)/dx);
        int y_idx_left = (int)((y - d_max - Ymin)/dy);
        int y_idx_right = (int)((y + d_max - Ymin)/dy);
        x_idx_left = max(0, x_idx_left);
        x_idx_right = min(Maxx, x_idx_right);
        y_idx_left = max(0, y_idx_left);
        y_idx_right = min(Maxy, y_idx_right);
        for (int ix = x_idx_left; ix < x_idx_right; ix++) {
            double xg = Xmin + ix*dx;
            for (int iy = y_idx_left; iy < y_idx_right; iy++) {
                double yg = Ymin + iy*dy;
                double dc = (x-xg)*(x-xg) + (y-yg)*(y-yg);
                if (shape_of_nucleons == 1) {
                    // "Checker" nucleons:
                    if (dc>dsq) continue;
                    if (participant[ipart]->isNucl() == 1) {
                        TA1[ix][iy] += areai;
                    } else if (participant[ipart]->isNucl() == 2) {
                        TA2[ix][iy] += areai;
                    } else {
                        cout << " Error in getTA2() " << endl;
                        exit(1);
                    }
                } else if (shape_of_nucleons >= 2
                           && shape_of_nucleons <= 9) {
                    // Gaussian nucleons:
                    // skip distant nucleons, speeds things up; 
                    // one may need to relax
                    if (dc > dc_sq_max_gaussian) continue;
                    double density = GaussianNucleonsCal::get2DHeightFromWidth(nucleon_width)*exp(-dc/(2.*nucleon_width*nucleon_width)); // width given from GaussianNucleonsCal class, height from the requirement that density should normalized to 1
                    if (participant[ipart]->isNucl() == 1) {
                        TA1[ix][iy] += density;
                    } else if (participant[ipart]->isNucl() == 2) {
                        TA2[ix][iy] += density;
                    } else {
                        cout << " Error in getTA2() " << endl;
                        exit(1);
                    }
                }
            }
        }
    }
    
    for (int ix = 0; ix < Maxx; ix++) {
        for (int iy = 0; iy < Maxy; iy++) {
            // keep track of highest density
            int i = int((TA1[ix][iy])/areai+0.5);
            int j = int((TA2[ix][iy])/areai+0.5);
            imax = max(imax,i);
            imax = max(imax,j);
        }
    }
    if (which_mc_model == 1) {
        if (imax>=tmax) {
            cout << "# WARNING: in MCnucl::getTA2() : imax=" << imax
                 << " should be less than tmax=" << tmax << endl;
            exit(0);
        }
    }
    delete [] part_x;
    delete [] part_y;
}

// calculate the binary collision density in the transverse plane
void MCnucl::calculate_rho_binary()
{
   int ncoll=binaryCollision.size();
   double dc_sq_max_gaussian = 25.*entropy_gaussian_width_sq;
   double d_max;
   if (shape_of_nucleons == 1)
       d_max = 2.*sqrt(dsq);
   if (shape_of_nucleons >= 2 && shape_of_nucleons <=9)
       d_max = 5.*entropy_gaussian_width;
   double *binary_x, *binary_y;
   binary_x = new double [ncoll];
   binary_y = new double [ncoll];
   for(int icoll = 0; icoll < ncoll; icoll++)
   {
       binary_x[icoll] = binaryCollision[icoll]->getX();
       binary_y[icoll] = binaryCollision[icoll]->getY();
   }
   for(int ir = 0; ir < Maxx; ir++)
       for(int jr = 0; jr < Maxy; jr++)
           rho_binary[ir][jr] = 0.0;

   for(int icoll = 0; icoll < ncoll; icoll++)
   {
       double x = binary_x[icoll];
       double y = binary_y[icoll];
       int x_idx_left = (int)((x - d_max - Xmin)/dx);
       int x_idx_right = (int)((x + d_max - Xmin)/dx);
       int y_idx_left = (int)((y - d_max - Ymin)/dy);
       int y_idx_right = (int)((y + d_max - Ymin)/dy);
       x_idx_left = max(0, x_idx_left);
       x_idx_right = min(Maxx, x_idx_right);
       y_idx_left = max(0, y_idx_left);
       y_idx_right = min(Maxy, y_idx_right);
       for(int ir = x_idx_left; ir < x_idx_right; ir++)
       {
          double xg = Xmin + ir*dx;
          for(int jr = y_idx_left; jr < y_idx_right; jr++)
          {
             double yg = Ymin + jr*dy;
             double dc = (x-xg)*(x-xg) + (y-yg)*(y-yg);
             if (shape_of_nucleons == 1)
             {
                if(dc <= dsq) 
                   rho_binary[ir][jr] += (10.0/siginNN);
             }
             else if (shape_of_nucleons>=2 && shape_of_nucleons<=9)
             {
                if (dc > dc_sq_max_gaussian) 
                   continue; // skip small numbers to speed up
                rho_binary[ir][jr] += GaussianNucleonsCal::get2DHeightFromWidth(entropy_gaussian_width)*exp(-dc/(2*entropy_gaussian_width_sq)); 
                // this density is normalized to 1, to be consistent with the disk-like treatment; 
             }
          }
       }
   }
   delete [] binary_x;
   delete [] binary_y;
}

// calculate the spectator density in the transverse plane
void MCnucl::calculate_spectator_density()
{
   int n_nucleon = spectators.size();
   double dc_sq_max_gaussian = 25.*entropy_gaussian_width_sq;
   double d_max;
   if (shape_of_nucleons == 1)
       d_max = 2.*sqrt(dsq);
   if (shape_of_nucleons >= 2 && shape_of_nucleons <=9)
       d_max = 5.*entropy_gaussian_width;
   double *spectator_x, *spectator_y, *spectator_rap;
   spectator_x = new double [n_nucleon];
   spectator_y = new double [n_nucleon];
   spectator_rap = new double [n_nucleon];
   for(int inucleon = 0; inucleon < n_nucleon; inucleon++)
   {
       spectator_x[inucleon] = spectators[inucleon]->getX();
       spectator_y[inucleon] = spectators[inucleon]->getY();
       spectator_rap[inucleon] = spectators[inucleon]->getRapidity_Y();
   }
   for(int ir = 0; ir < Maxx; ir++)
   {
       for(int jr = 0; jr < Maxy; jr++)
       {
           spectator_1[ir][jr] = 0.0;
           spectator_2[ir][jr] = 0.0;
       }
   }

   for(int inucleon = 0; inucleon < n_nucleon; inucleon++)
   {
       double x = spectator_x[inucleon];
       double y = spectator_y[inucleon];
       int nucleus_id = 0;
       if(spectator_rap[inucleon] > 0.)
           nucleus_id = 1;
       else
           nucleus_id = 2;

       int x_idx_left = (int)((x - d_max - Xmin)/dx);
       int x_idx_right = (int)((x + d_max - Xmin)/dx);
       int y_idx_left = (int)((y - d_max - Ymin)/dy);
       int y_idx_right = (int)((y + d_max - Ymin)/dy);
       x_idx_left = max(0, x_idx_left);
       x_idx_right = min(Maxx, x_idx_right);
       y_idx_left = max(0, y_idx_left);
       y_idx_right = min(Maxy, y_idx_right);
       for(int ir = x_idx_left; ir < x_idx_right; ir++)
       {
          double xg = Xmin + ir*dx;
          for(int jr = y_idx_left; jr < y_idx_right; jr++)
          {
             double yg = Ymin + jr*dy;
             double dc = (x-xg)*(x-xg) + (y-yg)*(y-yg);
             if (shape_of_nucleons == 1)
             {
                if(dc <= dsq) 
                {
                   if(nucleus_id == 1)
                       spectator_1[ir][jr] += (10.0/siginNN);
                   else
                       spectator_2[ir][jr] += (10.0/siginNN);
                }
             }
             else if (shape_of_nucleons>=2 && shape_of_nucleons<=9)
             {
                if (dc > dc_sq_max_gaussian) 
                   continue; // skip small numbers to speed up

                if(nucleus_id == 1)
                    spectator_1[ir][jr] += GaussianNucleonsCal::get2DHeightFromWidth(entropy_gaussian_width)*exp(-dc/(2*entropy_gaussian_width_sq)); 
                    // this density is normalized to 1, to be consistent with the disk-like treatment; 
                 else
                    spectator_2[ir][jr] += GaussianNucleonsCal::get2DHeightFromWidth(entropy_gaussian_width)*exp(-dc/(2*entropy_gaussian_width_sq)); 
             }
          }
       }
   }
   delete [] spectator_x;
   delete [] spectator_y;
   delete [] spectator_rap;
}

double MCnucl::get_spectator_density(int nucleus_id, int x, int y)
{
    if(nucleus_id == 1)
        return(spectator_1[x][y]);
    else
        return(spectator_2[x][y]);
}

// --- initializes dN/dyd2rt (or dEt/...) on 2d grid for rapidity slice iy
//     and integrates it to obtain dN/dy (or dEt/dy) ---
void MCnucl::setDensity(int iy, int ipt)
{
  // which_mc_model==1 -> KLN-like
  if (which_mc_model==1 && ipt>=0 && (dndydptTable==0)) {
    cout << "ERROR in MCnucl::setDensity() : pt-bin=" << ipt <<
      " but no dndydptTable !" << endl;
    exit(0);
  }

  // which_mc_model==1 -> KLN-like
  if (which_mc_model==1 && ipt<0 && (dndyTable==0)) {
    cout <<
     "ERROR in MCnucl::setDensity() : pt-integrated yields require dndyTable !" << endl;
    exit(0);
  }

  double tblmax=0, table_result=0;

  rapidity=rapMin + (rapMax-rapMin)/binRapidity*iy;
  dndy=0.0;
 
  double dc_sq_max_gaussian = 25.*entropy_gaussian_width_sq;
  double d_max;
  if (shape_of_nucleons == 1)
      d_max = 2.*sqrt(dsq);
  if (shape_of_nucleons >= 2 && shape_of_nucleons <=9)
      d_max = 5.*entropy_gaussian_width;
  int npart = participant.size();
  double *part_x = new double [npart];
  double *part_y = new double [npart];
  for(int i = 0; i < npart; i++)
  {
      part_x[i] = participant[i]->getX();
      part_y[i] = participant[i]->getY();
  }
  int ncoll = binaryCollision.size();
  double *binary_x = new double [ncoll];
  double *binary_y = new double [ncoll];
  for(int i = 0; i < ncoll; i++)
  {
      binary_x[i] = binaryCollision[i]->getX();
      binary_y[i] = binaryCollision[i]->getY();
  }

  if(which_mc_model==1) // MC-KLN
  {
    for(int ir=0;ir<Maxx;ir++)  // loop over 2d transv. grid
      for(int jr=0;jr<Maxy;jr++) 
      {
        double xg = Xmin + ir*dx;
        double yg = Ymin + jr*dy;
        double di = TA1[ir][jr]/dT;
        double dj = TA2[ir][jr]/dT;
        tblmax=max(di,tblmax);
        tblmax=max(dj,tblmax);
        if((di<0 || di>=tmax-2) || (dj<0 || dj>=tmax-2) ) {
            cerr << "di= " << di << " dj= " << dj
                << " You should increase the dimension of dndyTable"
                << endl;
            exit(1);
        }
        int i = floor(di); int j = floor(dj);
        if (ipt<0) // without pt dependence
        {
          table_result = sixPoint2dInterp(di-i, dj-j, // x and y value, in lattice unit (dndyTable_step -> 1)
          dndyTable[iy][i][j], dndyTable[iy][i][j+1], dndyTable[iy][i][j+2], dndyTable[iy][i+1][j], dndyTable[iy][i+1][j+1], dndyTable[iy][i+2][j]);
          rho->setDensity(iy,ir,jr,table_result);
        }
        else // with pt dependence
        {
          table_result = sixPoint2dInterp(di-i, dj-j, // x and y value, in lattice unit (dndyTable_step -> 1)
          dndydptTable[iy][i][j][ipt], dndydptTable[iy][i][j+1][ipt], dndydptTable[iy][i][j+2][ipt], dndydptTable[iy][i+1][j][ipt], dndydptTable[iy][i+1][j+1][ipt], dndydptTable[iy][i+2][j][ipt]);
          rho->setDensity(iy,ir,jr, ipt, table_result);
        }
        //dndy += dndyTable[iy][i][j];
        dndy += table_result;
      }
  }
  else if(which_mc_model==5) // MC-Glb, entropy centered at wounded nucleons and binary collision centers
  {
      double **rhop, **tab;
      rhop = new double * [Maxx];
      tab = new double * [Maxx];
      for(int ir = 0; ir < Maxx; ir++)
      {
          rhop[ir] = new double [Maxy];
          tab[ir] = new double [Maxy];
          for(int jr = 0; jr < Maxy; jr++)
          {
              rhop[ir][jr] = 0.0;
              tab[ir][jr] = 0.0;
              sd_TA1[ir][jr] = 0.0;
              sd_TA2[ir][jr] = 0.0;
          }
      }
      // wounded nucleon treatment:
      if (sub_model==1) // "classical" Glb
      {
          //rhop = (TA1[ir][jr]+TA2[ir][jr])*(1.0-Alpha)/2;
          double fluctfactor = 1.0;
          for(unsigned int ipart=0; ipart<participant.size(); ipart++) 
          {
            double x = part_x[ipart];
            double y = part_y[ipart];
            if (CCFluctuationModel > 5)
               fluctfactor = participant[ipart]->getfluctfactor();
            int x_idx_left = (int)((x - d_max - Xmin)/dx);
            int x_idx_right = (int)((x + d_max - Xmin)/dx);
            int y_idx_left = (int)((y - d_max - Ymin)/dy);
            int y_idx_right = (int)((y + d_max - Ymin)/dy);
            x_idx_left = max(0, x_idx_left);
            x_idx_right = min(Maxx, x_idx_right);
            y_idx_left = max(0, y_idx_left);
            y_idx_right = min(Maxy, y_idx_right);
            for(int ir = x_idx_left; ir < x_idx_right; ir++)
            {
               double xg = Xmin + ir*dx;
               for(int jr = y_idx_left; jr < y_idx_right; jr++)
               {
                   double yg = Ymin + jr*dy;
                   double dc = (x-xg)*(x-xg) + (y-yg)*(y-yg);
                   if (shape_of_entropy == 1) // "Checker" nucleons:
                   {
                     if(dc>dsq) continue;
                     double areai = 10.0/siginNN;
                     rhop[ir][jr] += fluctfactor*areai;
                     if (participant[ipart]->isNucl() == 1) {
                         sd_TA1[ir][jr] += fluctfactor*areai;
                     } else if (participant[ipart]->isNucl() == 2) {
                         sd_TA2[ir][jr] += fluctfactor*areai;
                     }
                   }
                   else if (shape_of_entropy>=2 && shape_of_entropy<=9) // Gaussian nucleons:
                   {
                     // skip distant nucleons, speeds things up; one may need to relax
                     if (dc>dc_sq_max_gaussian) continue;
                     double density;
                     if (shape_of_entropy == 3) {
                        density = participant[ipart]->getParticle()->getInternalStructDensity(xg, yg, quark_width, gaussDist); // width given from GaussianNucleonsCal class, height from the requirement that density should normalized to 1
                     } else {
                        density = fluctfactor*GaussianNucleonsCal::get2DHeightFromWidth(entropy_gaussian_width)*exp(-dc/(2.*entropy_gaussian_width_sq)); // width given from GaussianNucleonsCal class, height from the requirement that density should normalized to 1
                     }
                     rhop[ir][jr] += density;
                     if (participant[ipart]->isNucl() == 1) {
                         sd_TA1[ir][jr] += density;
                     } else if (participant[ipart]->isNucl() == 2) {
                         sd_TA2[ir][jr] += density;
                     }
                   }
               }
            }
          }
          double prefactor = (1.0 - Alpha)/2.;
          for(int ir = 0; ir < Maxx; ir++)
              for(int jr = 0; jr < Maxy; jr++)
                  rhop[ir][jr] = rhop[ir][jr]*prefactor;
      }
      else if (sub_model==2) // "Ulrich" Glb
      {
          for(int ir = 0; ir < Maxx; ir++)
              for(int jr = 0; jr < Maxy; jr++)
                  rhop[ir][jr] = 0.0; // no "wounded" contribution
      }
      else
      {
          cout << "MCnucl::setDensity error: which_mc_model is set to " << which_mc_model << ", but the associated sub_model " << sub_model << " is not recognized." << endl;
          exit(-1);
      }

      // binary collision treatment:
      if(Alpha > 1e-8)
      {
          int ncoll=binaryCollision.size();
          double fluctfactor = 1.0;
          for(int icoll=0;icoll<ncoll;icoll++)
          {
              double x = binary_x[icoll];
              double y = binary_y[icoll];
              if(CCFluctuationModel > 5)
                  fluctfactor = binaryCollision[icoll]->getfluctfactor();
              int x_idx_left = (int)((x - d_max - Xmin)/dx);
              int x_idx_right = (int)((x + d_max - Xmin)/dx);
              int y_idx_left = (int)((y - d_max - Ymin)/dy);
              int y_idx_right = (int)((y + d_max - Ymin)/dy);
              x_idx_left = max(0, x_idx_left);
              x_idx_right = min(Maxx, x_idx_right);
              y_idx_left = max(0, y_idx_left);
              y_idx_right = min(Maxy, y_idx_right);
              for(int ir = x_idx_left; ir < x_idx_right; ir++)
              {
                 double xg = Xmin + ir*dx;
                 for(int jr = y_idx_left; jr < y_idx_right; jr++)
                 {
                     double yg = Ymin + jr*dy;
                     double dc = (x-xg)*(x-xg) + (y-yg)*(y-yg);
                     if (shape_of_entropy==1)
                     {
                         if(dc <= dsq) 
                             tab[ir][jr] += fluctfactor*(10.0/siginNN)*(Alpha + (1.-Alpha)*binaryCollision[icoll]->additional_weight); // second part in the paranthesis is for Uli-Glb model
                     }
                     else if (shape_of_entropy>=2 && shape_of_entropy<=9)
                     {
                         if (dc <= dc_sq_max_gaussian)
                             tab[ir][jr] += fluctfactor*GaussianNucleonsCal::get2DHeightFromWidth(entropy_gaussian_width)*exp(-dc/(2*entropy_gaussian_width_sq))*(Alpha + (1.-Alpha)*binaryCollision[icoll]->additional_weight); // this density is normalized to 1, to be consisitent with the disk-like treatment; second part in the last parathesis is for Uli-Glb model
                     }
                 }
              }
          }
      } 
      else
      {
          for(int ir = 0; ir < Maxx; ir++)
              for(int jr = 0; jr < Maxy; jr++)
                  tab[ir][jr] = 0.;
      } // if(Alpha > 1e-8)

      // set entropy density
      for(int ir = 0; ir < Maxx; ir++)
          for(int jr = 0; jr < Maxy; jr++)
          {
              double density = rhop[ir][jr] + tab[ir][jr]; // "rhop" (n_WN) and "tab" (n_bin) have already been multiplied by (1-x)/2 and x (x=Alpha here) correspondingly
              rho->setDensity(iy,ir,jr,density);
              dndy += density;
          }
      for(int ir = 0; ir < Maxx; ir++)
      {
          delete [] rhop[ir];
          delete [] tab[ir];
      }
      delete [] rhop;
      delete [] tab;
  }
  if(dndy<1e-15)
  {
    cout << "MCnucl::setDensity dndy = 0 !!  y= " << rapidity
     << " dndy= " << dndy << endl;
    exit(0);
  }

  // Should I include additional fluctuation for MCKLN?
  if (CCFluctuationModel>0 && CCFluctuationModel <= 5) fluctuateCurrentDensity(iy);

  delete [] part_x;
  delete [] part_y;
  delete [] binary_x;
  delete [] binary_y;
}


//----------------------------------------------------------------------
void MCnucl::fluctuateCurrentDensity(int iy)
// Fluctuate the density profile in rho
{
    if (CCFluctuationModel==1) // use constant k/n:
    {
        for(int ir=0;ir<Maxx;ir++)  // loop over 2d transv. grid
        for(int jr=0;jr<Maxy;jr++)
        {
            double nb = rho->getDensity(iy,ir,jr)*dx*dy;
            double n = nbd->rand(nb/(nb+CCFluctuationK), CCFluctuationK);
            rho -> setDensity(iy,ir,jr,n/(dx*dy));
        }
    }
    else if (CCFluctuationModel==2)
    {
        // double kpp = 1.0/M_PI*dx*dy*1.0*kln->getLambdaQCD(); // second 1.0 is delta-eta
        double kpp = 1.0/M_PI*dx*dy*1.0*(0.25*0.25/HBARC/HBARC); // lambdaQCD=0.25/hbarC
        for(int ir=0;ir<Maxx;ir++)  // loop over 2d transv. grid
        for(int jr=0;jr<Maxy;jr++)
        {
            double k = kpp*min(TA1[ir][jr], TA2[ir][jr])*siginNN/10;
            double nb = rho->getDensity(iy,ir,jr)*dx*dy;
            double n;
            if (nb<1e-10)
                n = nb;
            else
                n = nbd->rand(nb/(nb+k), k);
            rho -> setDensity(iy,ir,jr,n/(dx*dy));
        }
    }
    else
    {
        cout << "MCnucl::fluctuateCurrentDensity error: CCFluctuationModel not supported." << endl;
        cout << "MCnucl:: CCFluctuationModel = " << CCFluctuationModel << endl;
        exit(-1);
    }
}


/***
   generates lookup table for dN/dy for varying proj/target thicknesses
***/
void MCnucl::makeTable()
{
  cout << "MCnucl::makeTable(): precalculating dNdy for all combinations of Ta and Tb." << endl;

  dT=10.0/siginNN/overSample;   // elementary thickness step
  if (shape_of_nucleons>=2 && shape_of_nucleons<=9) {  // Gaussian nucleons require finer steps
    tmax = paraRdr->getVal("tmax_subdivision")*(tmax-1) + 1;
    dT /= paraRdr->getVal("tmax_subdivision");
  }

  // allocate table
  dndyTable = new double** [binRapidity];
  for(int iy=0;iy<binRapidity;iy++) {
    dndyTable[iy] = new double* [tmax];
    for(int j=0;j<tmax;j++) dndyTable[iy][j] = new double [tmax];
  }

int progress_counter = 0, progress_percent = 0, last_update = 0;
//===========================================================================
  for(int iy=0;iy<binRapidity;iy++) { // loop over rapidity bins
    double y = rapMin+(rapMax-rapMin)/binRapidity*iy;

    for(int i=0;i<tmax;i++) {   // loop over proj thickness
      double ta1 = dT*i;
      for(int j=0;j<tmax;j++) { // loop over targ thickness
        double ta2 = dT*j;
        if(i>0 && j>0) {  // store corresponding dN/dy in lookup table
          // small-x gluons via kt-factorization
          dndyTable[iy][i][j] = kln->getdNdy(y,ta1,ta2, -1, PT_order); 
          // add large-x partons via DHJ formula if required
          if (val)
            dndyTable[iy][i][j] += val->getdNdy(y,ta1,ta2);
          //cout << ta1 << ", " << ta2 << ", " << dndyTable[iy][i][j] << endl;
        } else dndyTable[iy][i][j] = 0.0;
      progress_counter++;
      progress_percent = (progress_counter*100) / (binRapidity*tmax*tmax);
      if(((progress_percent%10) == 0) && (progress_percent != last_update))
    {
       cout << progress_percent << "% : " << std::flush;
       last_update = progress_percent;
      }
      }
    }
  }
  cout << endl;
  //===========================================================================
  cout << "MCnucl::makeTable(): done" << endl << endl;

  dumpdNdyTable4Col("data/dNdyTable.dat", dndyTable, 0);
}



/***
   generates lookup table for dN/dyd2pt for varying proj/target thicknesses
***/
void MCnucl::makeTable(double ptmin, double dpt, int iPtmax)
{
  cout << "MCnucl::makeTable(double, double, int): generating dN/dyd2pt lookup table... " << endl;

  dT=10.0/siginNN/overSample;   // elementary thickness step
  if (shape_of_nucleons>=2 && shape_of_nucleons<=9) {  // Gaussian nucleons require finer steps
    tmax = paraRdr->getVal("tmax_subdivition")*(tmax -1 ) + 1;
    dT /= paraRdr->getVal("tmax_subdivition");
  }

  // range of thicknesses for pt dependent lookup table
  tmaxPt = tmax;
  //tmaxPt = tmax/5;   // sufficient for pp

  // allocate table
  iptmax = iPtmax;  // destructor needs to delete table, store size
  dndydptTable = new double*** [binRapidity];
  for (int iy=0; iy<binRapidity; iy++) {
    dndydptTable[iy] = new double** [tmaxPt];
    for (int j=0; j<tmaxPt; j++) {
      dndydptTable[iy][j] = new double* [tmaxPt];
      for (int i=0; i<tmaxPt; i++)
        dndydptTable[iy][j][i] = new double [iptmax];
    }
  }

  
  int progress_counter = 0, progress_percent = 0, last_update = 0;
  //===========================================================================

  for(int iy=0;iy<binRapidity;iy++) { // loop over rapidity bins
    double y = rapMin+(rapMax-rapMin)/binRapidity*iy;
    for (int i=0; i<tmaxPt; i++) {   // loop over proj thickness
      double ta1 = dT*i;
      for (int j=0; j<tmaxPt; j++) { // loop over targ thickness
        double ta2 = dT*j;
        if (val) val->lgXClearPtTable(); // recomp. large-rap pt-distributions
        for (int ipt=0; ipt<iptmax; ipt++) { // loop over pt bins
          if(i>0 && j>0) {  // store corresponding dN/dyd2pt in lookup table
            // high-rap *hadrons* via DHJ formula
            if (val)  dndydptTable[iy][i][j][ipt] =
            val->getdNdyd2pt(y,ta1,ta2,ptmin+ipt*dpt);
            else
              // small-x gluons via kt-factorization;  fixed pt, no integration
              dndydptTable[iy][i][j][ipt] = kln->getdNdy(y,ta1,ta2,ptmin+ipt*dpt);
          } else dndydptTable[iy][i][j][ipt] = 0.0;        
          progress_counter++;
          progress_percent = (progress_counter*100) / (binRapidity*tmax*tmax*iptmax);
          if(((progress_percent%10) == 0) && (progress_percent != last_update))
          {
           cout << progress_percent << "% : " << std::flush;
           last_update = progress_percent;
          }
        }
      }
    }
  }
  cout << "MCnucl::makeTable(double, double, int): done" << endl;
}

void MCnucl::dumpdNdyTable4Col(char filename[], double *** dNdyTable, const int iy)
{
  ofstream of;
  of.open(filename, std::ios_base::app);

  double y = rapMin+(rapMax-rapMin)/binRapidity*iy;

    for(int i=0;i<tmax;i++) {   // loop over proj thickness
      double ta1 = dT*i;
      for(int j=0;j<tmax;j++) { // loop over targ thickness
        double ta2 = dT*j;
        if(i>0 && j>0) {
          of << fixed << setprecision(3) << setw(10) <<  y
             << fixed << setprecision(3) << setw(10) <<  ta1
             << fixed << setprecision(3) << setw(10) <<  ta2
             << setprecision(12) << setw(22) <<  dNdyTable[iy][i][j]
             << endl;
        }
      }
    }
     of.close();
}

void MCnucl::dumpdNdydptTable5Col(char filename[], double **** dNdydptTable, const int iy)
{
  ofstream of;
  of.open(filename, std::ios_base::out);

  double y = rapMin+(rapMax-rapMin)/binRapidity*iy;
  double iptmax = MaxPT;

    for(int i=0;i<tmax;i++) {   // loop over proj thickness
      double ta1 = dT*i;
      for(int j=0;j<tmax;j++) { // loop over targ thickness
        double ta2 = dT*j;
        if(i>0 && j>0) {
          for(int ipt=0;ipt<iptmax; ipt++) { //loop over Pt
            double ptstep = PTmin+ipt*dpt;
          of << fixed << setprecision(3) << setw(10) <<  y
             << fixed << setprecision(3) << setw(10) <<  ta1
             << fixed << setprecision(3) << setw(10) <<  ta2
             << fixed << setprecision(3) << setw(10) <<  ptstep;
          of << scientific << setprecision(12) << setw(22) <<  dNdydptTable[iy][i][j][ipt]
             << endl;
             }
        }
      }
    }
     of.close();
}


void MCnucl::deleteNucleus()
{
  for(int i=0;i<(int)nucl1.size();i++) {
    delete nucl1[i];
  }
  for(int i=0;i<(int)nucl2.size();i++) {
    delete nucl2[i];
  }
  for(int i=0;i<(int)participant.size();i++) {
    delete participant[i];
  }
  for(int i=0;i<(int)spectators.size();i++) {
    delete spectators[i];
  }
  for(int i=0;i<(int)binaryCollision.size();i++) {
    delete binaryCollision[i];
  }
  nucl1.clear();
  nucl2.clear();
  participant.clear();
  spectators.clear();
  binaryCollision.clear();
}


double MCnucl::Angle(const double x,const double y)
{
    double angl=0.0;
    double r=sqrt(x*x+y*y);
    if(r < 1e-20) return angl;

    if(abs(x)/r < 0.8) {
        //angl=sign(acos(x/r),y)
        angl = y>0 ? abs(acos(x/r)): -abs(acos(x/r));
    }else {
        angl=asin(y/r);
        if(x < 0.0 && angl >= 0.0)
          angl=M_PI-angl;
        else if(x < 0.0)
          angl=-M_PI-angl;
    }

    return angl;
}

void MCnucl::recenterGrid(int iy, int n)
{
  rho->getCMAngle(iy, n);
  rho->recenterParticle(participant, binaryCollision, spectators, iy);
}

void MCnucl::rotateGrid(int iy, int n)
{
  rho->getCMAngle(iy, n);
  rho->recenterParticle(participant, binaryCollision, spectators, iy);
  rho->rotateParticle(participant, binaryCollision, spectators, iy);
}


void MCnucl::dumpBinaryTable(char filename[])
{
  double x,y;
  ofstream of;

  of.open(filename, std::ios_base::app);
  for (int idx=0; idx<binaryCollision.size(); idx++)
  {
    x = binaryCollision[idx]->getX();
    y = binaryCollision[idx]->getY();
    of  << setprecision(3) << setw(10) << x
        << setprecision(3) << setw(10) << y
        << endl;
  }
  of.close();
  
  /* for debug
  of.open("data/wounded.data");
  for (int idx=0; idx<participant.size(); idx++)
  {
    x = participant[idx]->getX();
    y = participant[idx]->getY();
    of  << setprecision(3) << setw(10) << x
        << setprecision(3) << setw(10) << y
        << endl;
  }
  of.close();

  of.open("data/nucl1.data");
  for (int idx=0; idx<nucl1.size(); idx++)
  {
    x = nucl1[idx]->getX();
    y = nucl1[idx]->getY();
    of  << setprecision(3) << setw(10) << x
        << setprecision(3) << setw(10) << y
        << endl;
  }
  of.close();

  of.open("data/nucl2.data");
  for (int idx=0; idx<nucl2.size(); idx++)
  {
    x = nucl2[idx]->getX();
    y = nucl2[idx]->getY();
    of  << setprecision(3) << setw(10) << x
        << setprecision(3) << setw(10) << y
        << endl;
  }
  of.close();
  */
}

void MCnucl::dumpparticipantTable(char filename[])
{
  double x,y;
  ofstream of;

  of.open(filename, std::ios_base::app);
  for (int idx=0; idx<participant.size(); idx++)
  {
    x = participant[idx]->getX();
    y = participant[idx]->getY();
    int id = participant[idx]->isNucl();
    of  << setprecision(3) << setw(10) << x << "   "
        << setprecision(3) << setw(10) << y << "   " 
        << id << endl;
  }
  of.close();
}


int MCnucl::getSpectators()
{
  int count = 0;
  double ecm = paraRdr->getVal("ecm");
  //calculate the rapidity_Y for spectators at a given collision energy
  double v_z = sqrt(1. - 1./((ecm/2.)*(ecm/2.)));
  double rapidity_Y = 0.5*log((1. + v_z)/(1. - v_z + 1e-100));
  for(int i=0;i<(int)nucl1.size();i++) { // loop over proj. nucleons
    if(nucl1[i]->getNumberOfCollision()==0) {
      double x1 = nucl1[i]->getX();
      double y1 = nucl1[i]->getY();
      spectators.push_back(new Spectator(x1, y1, rapidity_Y));
      count++;
    }
  }
  for(int i=0;i<(int)nucl2.size();i++) { // loop over targ. nucleons
    if(nucl2[i]->getNumberOfCollision()==0) {
      double x2 = nucl2[i]->getX();
      double y2 = nucl2[i]->getY();
      spectators.push_back(new Spectator(x2, y2, -rapidity_Y));
      count++;
    }
  }
  return(count);
}

void MCnucl::dumpSpectatorsTable(int event)
{
  double x, y, rap;
  ostringstream of_stream; 
  ofstream of;
  of_stream << "data/Spectators_event_" << event << ".dat";
  
  of.open(of_stream.str().c_str());
  for (int idx=0; idx<spectators.size(); idx++)
  {
    x = spectators[idx]->getX();
    y = spectators[idx]->getY();
    rap = spectators[idx]->getRapidity_Y();
    of  << scientific << setprecision(4) << setw(10) 
        << x << "  " << y << "  " << rap
        << endl;
  }
  of.close();
}

double MCnucl::sampleFluctuationFactorforParticipant()
{
   double eps = 1e-8;
   double fluctfactor = 1.0;
   double Gamma_k = 1./ccFluctuationGammaTheta;
   double k_part = (1 - Alpha + eps)/2.*Gamma_k;
   double theta_part = 2./(1 - Alpha + eps)*ccFluctuationGammaTheta;
   if(CCFluctuationModel == 6)  //Gamma distribution for MC-Glauber
      fluctfactor = gsl_ran_gamma(gslRng, k_part, theta_part);
   
   return(fluctfactor);
}

double MCnucl::sampleFluctuationFactorforBinaryCollision()
{
   double eps = 1e-8;
   double fluctfactor = 1.0;
   double Gamma_k = 1./ccFluctuationGammaTheta;
   double k_binary = (Alpha+eps)*Gamma_k;
   double theta_binary = 1./(Alpha+eps)*ccFluctuationGammaTheta;
   if(CCFluctuationModel == 6)  //Gamma distribution for MC-Glauber
      fluctfactor = gsl_ran_gamma(gslRng, k_binary, theta_binary);
   
   return(fluctfactor);
}
