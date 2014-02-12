#include <iomanip>
#include <algorithm>
#include "MakeDensity.h"
#include "Regge96.h"
#include "KLNfunc.h"
#include "rcBKfunc.h"
#include "ParamDefs.h"
#include "MCnucl.h"
#include "EOS.h"

#include "ParameterReader.h"
#include "ArraySaver.h"
#include "Stopwatch.h"

using namespace std;

//***************************************************************************
// calculates dN/dyd2rt for AA events, each event is rotated in transv. plane
// to principal axis frame (angle determined by y=0 slice)

MakeDensity::MakeDensity(ParameterReader *paraRdr_in)
{
  paraRdr = paraRdr_in; // where all parameters are1

  eos.loadEOSFromFile("s95p-PCE/EOS_converted.dat","s95p-PCE/coeff.dat");

  int which_mc_model = paraRdr->getVal("which_mc_model");
  finalFactor = paraRdr->getVal("finalFactor");

  int lgX = paraRdr->getVal("lgX");

  bmin = paraRdr->getVal("bmin");
  bmax = paraRdr->getVal("bmax");
  Npmin = paraRdr->getVal("Npmin");
  Npmax = paraRdr->getVal("Npmax");

  cutdSdy = paraRdr->getVal("cutdSdy");
  cutdSdy_lowerBound = paraRdr->getVal("cutdSdy_lowerBound");
  cutdSdy_upperBound = paraRdr->getVal("cutdSdy_upperBound");

  // renaming of A for proj+targ
  Anucl1 = paraRdr->getVal("Aproj");
  Anucl2 = paraRdr->getVal("Atarg");

  // fix grid properties
  Xmax = paraRdr->getVal("maxx");
  Ymax = paraRdr->getVal("maxy");
  Xmin = -Xmax;
  Ymin = -Ymax;
  dx = paraRdr->getVal("dx");
  dy = paraRdr->getVal("dy");
  Maxx=(int)((Xmax-Xmin)/dx+0.1)+1;
  Maxy=(int)((Ymax-Ymin)/dy+0.1)+1;

  // rapidity binning
  binRapidity = paraRdr->getVal("ny");
  rapMin = -paraRdr->getVal("ymax");
  paraRdr->setVal("rapMin", rapMin);
  rapMax = -rapMin;
  paraRdr->setVal("rapMax", rapMax);

  // NN cross sections
  ecm = paraRdr->getVal("ecm");
  double sig = hadronxsec::totalXsection(ecm,0);
  double sigel = hadronxsec::elasticXsection(sig,ecm,0,0);
  siginNN = sig - sigel;

  // Unintegrated PT optns
  PTinte = paraRdr->getVal("PT_flag");
  PTmax  = paraRdr->getVal("PT_Max");
  PTmin  = paraRdr->getVal("PT_Min");
  dpt = paraRdr->getVal("d_PT");
  MaxPT=(int)((PTmax-PTmin)/dpt+0.1)+1;
  MixedMode = paraRdr->getVal("mixedMode"); 

  //operator option
  Operation = paraRdr->getVal("operation");
  // generate proj+targ MC thickness functions
  proj = new OverLap(Anucl1,siginNN,paraRdr->getVal("proj_deformed"));
  targ = new OverLap(Anucl2,siginNN,paraRdr->getVal("targ_deformed"));

  // overlap proj+targ on transverse grid
  mc = new MCnucl(paraRdr);

  // fix range of accepted Npart (centrality cut)
  mc->setCentralityCut((int) Npmin,(int) Npmax);

  // set KLN properties
  kln = 0;
  val = 0;
  wf = 0;
  if (which_mc_model==1) // MC-KLN
  {
    int mode = paraRdr->getVal("sub_model");
    if (mode >= rcBK) wf = new rcBKfunc(mode); // stupid coding style
    else wf = new KLNfunc();

    kln = new KLNModel(ecm, mode,wf);
    mc->setKLN(kln);
    // Get lambda value for kln
    kln->setLambdaSaturation(paraRdr->getVal("lambda"));

    // large-xF partons (rcBK mode only)
    if (lgX && mode >= rcBK) {
      val = new Large_x(kln,wf);  mc->setLgX(val);
    }

    // set overall normalization (0.735 for gaussian nucl, 0.7 for disk nucl.) this factor does not affect the initial condition the code outputs (i.e. gluon density)
    kln->setNormalization(1.0);

    // generates look-up table of dN/dy as fct of proj/targ thicknesses
    if(Operation==1 or Operation==3)   //generate event by event profile
    {
      if(PTinte>0 and MixedMode<=0)
        mc->makeTable();
      else if(PTinte<0)
        mc->makeTable(PTmin, dpt, MaxPT);  //Generate PT-unintegrated dNdydpt Table
      else if(MixedMode>0)   //output pT-integrated and unintegrated table at the same time
      {
        cout << "Mixed mode: Calculating two tables!" << endl;
        mc->makeTable();
        mc->makeTable(PTmin, dpt, MaxPT);  //Generate PT-unintegrated dNdydpt Table
      }
  
    }
    else if(Operation==9 or Operation==2)  //generate eccentricity table
    {
      if(PTinte>0)
        mc->makeTable();
      else
      {
        cout << "This option is not supported: " 
             << "Operation=" << Operation 
             << ", PTinte=" <<PTinte << endl;
        exit(0);
      }
    }
  }//if (which_mc_model==1)

}


MakeDensity::~MakeDensity()
{
  if (val) delete val;
  delete proj; delete targ;
  if (kln) delete kln;
  if (wf) delete wf;
  if (mc) delete mc;

}
//----------------------------------------------------------------------

//**********************************************************************
void MakeDensity::generate_profile_ebe_Jet(int nevent)
{
  // binary profile:
  char file_binary[] = "data/BinaryCollisionTable_event_%d.dat";

  // entropy profile:
  char file1_ecc[] = "data/sn_ecc_eccp_%%d_event_%d.dat";
  char file1_4col[] = "data/sd_event_%d_4col.dat";
  char file1_block[] = "data/sd_event_%d_block.dat";
  double *** dens1  = new double** [binRapidity];
  for(int iy=0;iy<binRapidity;iy++) {
    dens1[iy] =  new double* [Maxx]();
    for(int i=0;i<Maxx;i++) {
        dens1[iy][i] = new double[Maxy]();
        for (int j=0;j<Maxy;j++) dens1[iy][i][j]=0;
    }
  }

  // energy profile:
  char file2_ecc[] = "data/en_ecc_eccp_%%d_event_%d.dat";
  char file2_4col[] = "data/ed_event_%d_4col.dat";
  char file2_block[] = "data/ed_event_%d_block.dat"; 
  double *** dens2  = new double** [binRapidity];
  for(int iy=0;iy<binRapidity;iy++) {
    dens2[iy] =  new double* [Maxx]();
    for(int i=0;i<Maxx;i++) {
        dens2[iy][i] = new double[Maxy];
        for (int j=0;j<Maxy;j++) dens2[iy][i][j]=0;
    }
  }

  int from_order = paraRdr->getVal("ecc_from_order");
  int to_order = paraRdr->getVal("ecc_to_order");

  int use_sd = paraRdr->getVal("use_sd");
  int use_ed = paraRdr->getVal("use_ed");
  int use_block = paraRdr->getVal("use_block");
  int use_4col = paraRdr->getVal("use_4col");

  char buffer[200];
  double b;

  // event start.
  int event=1;
  while (event<=nevent)
  {
    int tries = 0;
    int binary = 0;

    while (binary==0 || mc->CentralityCut()==0)
    {
      b = sqrt((bmax*bmax - bmin*bmin)*drand48()+bmin*bmin);
      mc->generateNucleus(b,proj,targ);
      binary = mc->getBinaryCollision();
      //mc->dumpBinaryTable(); // for debugging
      if(binary==0 || mc->CentralityCut()==0) mc->deleteNucleus();
      tries++;
      if (tries>300)
      {
          cout << "===== MakeDensity::generate error =====" << endl
              << "No collisions detected after maximum number of tries." << endl
              << "What impact parameter are you using?" << endl;
          //exit(-1); // no more tries: something must be wrong.
          tries=0;
      }
    }
    Npart = mc->getNpart1()+mc->getNpart2();
    int Spectator = mc->getSpectators();
    mc->dumpSpectatorsTable(event);
    sprintf(buffer, file_binary, event);
    mc->dumpBinaryTable(buffer); // for collision profile

    // is event-by-event calculation; no need to rotate
    // compute density before rotation.
    mc->getTA2();

    bool cutdSdypassFlag = true;
    for(int iy=0;iy<binRapidity;iy++) {
      mc->setDensity(iy, -1);
      // cut total entropy
      if(iy == 0 && cutdSdy == 1)
      {
         double totaldSdy = gettotaldSdy(iy);
         if(totaldSdy < cutdSdy_lowerBound || totaldSdy > cutdSdy_upperBound)
         {
            cutdSdypassFlag = false;
            break;
         }
      }
      // output entropy profile
      if (use_sd)
      {
        setSd(dens1, iy); // includes factor multiplication
        sprintf(buffer, file1_ecc, event);
        dumpEccentricities(buffer, dens1, iy, from_order, to_order, Npart, binary, b);
        if (use_4col)
          {
            sprintf(buffer,file1_4col,event);
            dumpDensity4Col(buffer, dens1, iy);
          }
        if (use_block)
          {
            sprintf(buffer,file1_block,event);
            dumpDensityBlock(buffer, dens1, iy);
          }
      }
      // output energy profile
      if (use_ed)
      {
        setEd(dens2, iy); // includes factor multiplication
        sprintf(buffer,file2_ecc,event);
        dumpEccentricities(buffer, dens2, iy, from_order, to_order, Npart, binary, b);
        if (use_4col)
          {
            sprintf(buffer,file2_4col,event);
            dumpDensity4Col(buffer, dens2, iy);
          }
        if (use_block)
          {
            sprintf(buffer,file2_block,event);
            dumpDensityBlock(buffer, dens2, iy);
          }
      }
    } // <-> for(int iy=0;iy<binRapidity;iy++)

    mc->deleteNucleus();
    if(cutdSdypassFlag)
      event++;

  //break; # for debugging
  } // <-> while(event<events)

  // clean up
  for(int iy=0;iy<binRapidity;iy++) {
    for(int i=0;i<Maxx;i++) delete [] dens1[iy][i];
    delete [] dens1[iy];
  }
  delete [] dens1;

  for(int iy=0;iy<binRapidity;iy++) {
    for(int i=0;i<Maxx;i++) delete [] dens2[iy][i];
    delete [] dens2[iy];
  }
  delete [] dens2;
}
//----------------------------------------------------------------------

//**********************************************************************
void MakeDensity::generate_profile_ebe(int nevent)
{
  // entropy profile:
  char file1_4col[] = "data/sd_event_%d_4col.dat";
  char file1_block[] = "data/sd_event_%d_block.dat";
  char file1_5col[] = "data/sd_event_%d_5col.dat";
  char file1_ptCol[] = "data/sd_event_%d_ptCol.dat";    
  double *** dens1  = new double** [binRapidity];
  for(int iy=0;iy<binRapidity;iy++) {
    dens1[iy] =  new double* [Maxx]();
    for(int i=0;i<Maxx;i++) {
        dens1[iy][i] = new double[Maxy]();
        for (int j=0;j<Maxy;j++) dens1[iy][i][j]=0;
    }
  }

  double **** dens1pt  = new double*** [binRapidity];
    for(int iy=0;iy<binRapidity;iy++) {
    dens1pt[iy] =  new double** [Maxx]();
    for(int i=0;i<Maxx;i++) {
        dens1pt[iy][i] = new double* [Maxy]();
        for (int j=0;j<Maxy;j++) {
          dens1pt[iy][i][j]= new double [MaxPT]();
            for(int ipt=0;ipt<MaxPT;ipt++)
          dens1pt[iy][i][j][ipt]=0;}
    }
  }

  // energy profile:
  char file2_4col[] = "data/ed_event_%d_4col.dat";
  char file2_block[] = "data/ed_event_%d_block.dat";
  char file2_5col[] = "data/ed_event_%d_5col.dat";
  char file2_ptCol[] = "data/ed_event_%d_ptCol.dat";  
  double *** dens2  = new double** [binRapidity];
  for(int iy=0;iy<binRapidity;iy++) {
    dens2[iy] =  new double* [Maxx]();
    for(int i=0;i<Maxx;i++) {
        dens2[iy][i] = new double[Maxy];
        for (int j=0;j<Maxy;j++) dens2[iy][i][j]=0;
    }
  }
  double **** dens2pt  = new double*** [binRapidity];
    for(int iy=0;iy<binRapidity;iy++) {
    dens2pt[iy] =  new double** [Maxx]();
    for(int i=0;i<Maxx;i++) {
        dens2pt[iy][i] = new double* [Maxy]();
        for (int j=0;j<Maxy;j++) {
          dens2pt[iy][i][j]= new double [MaxPT]();
            for(int ipt=0;ipt<MaxPT;ipt++)
          dens2pt[iy][i][j][ipt]=0;}
    }
  }

  int use_sd = paraRdr->getVal("use_sd");
  int use_ed = paraRdr->getVal("use_ed");
  int use_block = paraRdr->getVal("use_block");
  int use_4col = paraRdr->getVal("use_4col");
  
  int use_5col = 0;
  int use_ptCol = 0;
  // sub-switch
  if(PTinte>0)
  {
    int use_5col = paraRdr->getVal("use_5col");
    int use_ptCol = paraRdr->getVal("use_ptCol");
  }

  char buffer[200];
  double b;

  // event start.
  int event=1;
  while (event<=nevent)
  {
    int tries = 0;
    int binary = 0;

    while (binary==0 || mc->CentralityCut()==0)
    {
      b = sqrt((bmax*bmax - bmin*bmin)*drand48()+bmin*bmin);
      mc->generateNucleus(b,proj,targ);
      binary = mc->getBinaryCollision();
      //mc->dumpBinaryTable(); // for debugging
      if(binary==0 || mc->CentralityCut()==0) mc->deleteNucleus();
      tries++;
      if (tries>300)
      {
          cout << "===== MakeDensity::generate error =====" << endl
              << "No collisions detected after maximum number of tries." << endl
              << "What impact parameter are you using?" << endl;
          //exit(-1); // no more tries: something must be wrong.
          tries=0;
      }
    }
    Npart = mc->getNpart1()+mc->getNpart2();
    //mc->dumpBinaryTable(); // for collision profile

    // is event-by-event calculation; no need to rotate
    // compute density before rotation.
    mc->getTA2();

    bool cutdSdypassFlag = true;
    if(PTinte>0)
    {
      for(int iy=0;iy<binRapidity;iy++) {
        mc->setDensity(iy, -1);
        // cut total entropy
        if(iy == 0 && cutdSdy == 1)
        {
           double totaldSdy = gettotaldSdy(iy);
           if(totaldSdy < cutdSdy_lowerBound || totaldSdy > cutdSdy_upperBound)
           {
              cutdSdypassFlag = false;
              break;
           }
        }
        // output entropy profile
        if (use_sd)
        {
          setSd(dens1, iy); // includes factor multiplication
          if (use_4col)
            {
              sprintf(buffer,file1_4col,event);
              dumpDensity4Col(buffer, dens1, iy);
            }
          if (use_block)
            {
              sprintf(buffer,file1_block,event);
              dumpDensityBlock(buffer, dens1, iy);
            }
        }
        // output energy profile
        if (use_ed)
        {
          setEd(dens2, iy); // includes factor multiplication
          if (use_4col)
            {
              sprintf(buffer,file2_4col,event);
              dumpDensity4Col(buffer, dens2, iy);
            }
          if (use_block)
            {
              sprintf(buffer,file2_block,event);
              dumpDensityBlock(buffer, dens2, iy);
            }
        }
      } // <-> for(int iy=0;iy<binRapidity;iy++)
      /* comment the following lines to let dNdyTable and dNdydPtTable use the same configuration*/   
      if(MixedMode<=0)   //keep the current configuration if need pT-unintegrated table
      {
        mc->deleteNucleus();
        if(cutdSdypassFlag)
          event++;
      }
  }// <-> if(PTinte>0) 

  if((PTinte<0 or MixedMode>0) and cutdSdypassFlag)
  //calculate pT unintegrated particle distribution
  {               
      for(int iy=0;iy<binRapidity;iy++) 
      {
        for(int ipt=0;ipt<MaxPT;ipt++) 
          mc->setDensity(iy,ipt);
      // output entropy profile
        if (use_sd)
        {
          setSd(dens1pt, iy); // includes factor multiplication
          if (use_5col)
          {
            sprintf(buffer,file1_5col,event);
            dumpDensity5Col(buffer, dens1pt, iy);
          }
          if (use_ptCol)     // jia test
          {
            sprintf(buffer,file1_ptCol,event);
            dumpDesityptCol(buffer, dens1pt, iy);
          }
        }
      
      // output energy profile
      if (use_ed)  //not applicable, since 
      {
        setEd(dens2pt, iy); // includes factor multiplication
        if (use_5col)
          {
            sprintf(buffer,file2_5col,event);
            dumpDensity5Col(buffer, dens2pt, iy);
          }
          if (use_ptCol)     // jia test
          {
            sprintf(buffer,file2_ptCol,event);
            dumpDesityptCol(buffer, dens2pt, iy);
          }
      }
     //<-> for for(int ipt=0;ipt<MaxPT;ipt++)
    } // <-> for(int iy=0;iy<binRapidity;iy++)
  
    mc->deleteNucleus();
    event++;
    } // <-> if(PTinte<0) 
  //break; # for debugging
  } // <-> while(event<events)

  // clean up
  for(int iy=0;iy<binRapidity;iy++) {
    for(int i=0;i<Maxx;i++) delete [] dens1[iy][i];
    delete [] dens1[iy];
  }
  delete [] dens1;

  for(int iy=0;iy<binRapidity;iy++) {
    for(int i=0;i<Maxx;i++) delete [] dens2[iy][i];
    delete [] dens2[iy];
  }
  delete [] dens2;

//clean up pt-dependence densities
  for(int iy=0;iy<binRapidity;iy++) {
    for(int i=0;i<Maxx;i++) {
      for(int j=0;j<Maxy;j++)
        delete [] dens1pt[iy][i][j];
      delete [] dens1pt[iy][i];
    }
    delete [] dens1pt[iy];
  }
  delete [] dens1pt;


  for(int iy=0;iy<binRapidity;iy++) {
    for(int i=0;i<Maxx;i++) {
      for(int j=0;j<Maxy;j++)
        delete [] dens2pt[iy][i][j];
      delete [] dens2pt[iy][i];
    }
    delete [] dens2pt[iy];
  }
  delete [] dens2pt;  
}
//----------------------------------------------------------------------



//**********************************************************************
void MakeDensity::generate_profile_average(int nevent)
{
  int average_from_order = paraRdr->getVal("average_from_order"), average_to_order = paraRdr->getVal("average_to_order");
  int number_of_orders = average_to_order - average_from_order + 1;

  int use_sd = paraRdr->getVal("use_sd");
  int use_ed = paraRdr->getVal("use_ed");
  int use_block = paraRdr->getVal("use_block");
  int use_4col = paraRdr->getVal("use_4col");
  int use_5col = 0;
  if(PTinte>0)
  {
    int use_5col = paraRdr->getVal("use_5col");
  }
  char file1_5col[] = "data/sdAvg_order_%d_5col.dat";

  // entropy profile:
  char file1_4col[] = "data/sdAvg_order_%d_4col.dat";
  char file1_block[] = "data/sdAvg_order_%d_block.dat";
  double **** dens1  = new double*** [number_of_orders];
  for(int iorder=0; iorder<number_of_orders; iorder++) // iorder starts from 0
  {
    dens1[iorder] = new double** [binRapidity];
    for(int iy=0;iy<binRapidity;iy++) {
      dens1[iorder][iy] =  new double* [Maxx]();
      for(int i=0;i<Maxx;i++) {
          dens1[iorder][iy][i] = new double[Maxy]();
          for (int j=0;j<Maxy;j++) dens1[iorder][iy][i][j]=0;
      }
    }
  }
  double ***** dens1pt = new double **** [number_of_orders];
  //entropy density for pt-unintegrated case: dens1pt(iorder, iy, x, y, ipt)
  for(int iorder=0; iorder<number_of_orders; iorder++) // iorder starts from 0
  {
    dens1pt[iorder] = new double*** [binRapidity];
    for(int iy=0;iy<binRapidity;iy++) {
      dens1pt[iorder][iy] =  new double** [Maxx]();
      for(int i=0;i<Maxx;i++) {
          dens1pt[iorder][iy][i] = new double* [Maxy]();
          for (int j=0;j<Maxy;j++) {
            dens1pt[iorder][iy][i][j] = new double[MaxPT]();
            for(int k=0;k<MaxPT;k++)
              dens1pt[iorder][iy][i][j][k] = 0.;
          }
      }
    }
  }

  // energy profile:
  char file2_4col[] = "data/edAvg_order_%d_4col.dat";
  char file2_block[] = "data/edAvg_order_%d_block.dat";
  double **** dens2  = new double*** [number_of_orders];
  char file2_5col[] = "data/edAvg_order_%d_5col.dat";
  for(int iorder=0; iorder<number_of_orders; iorder++) // iorder starts from 0
  {
    dens2[iorder] = new double** [binRapidity];
    for(int iy=0;iy<binRapidity;iy++) {
      dens2[iorder][iy] =  new double* [Maxx]();
      for(int i=0;i<Maxx;i++) {
          dens2[iorder][iy][i] = new double[Maxy]();
          for (int j=0;j<Maxy;j++) dens2[iorder][iy][i][j]=0;
      }
    }
  }
  //energy density for pt-unintegrated: dens2pt(iorder, iy, x, y, ipt)
  double *****dens2pt = new double**** [number_of_orders];
  for(int iorder=0; iorder<number_of_orders; iorder++) // iorder starts from 0
  {
    dens2pt[iorder] = new double*** [binRapidity];
    for(int iy=0;iy<binRapidity;iy++) {
      dens2pt[iorder][iy] =  new double** [Maxx]();
      for(int i=0;i<Maxx;i++) {
          dens2pt[iorder][iy][i] = new double* [Maxy]();
          for (int j=0;j<Maxy;j++) {
            dens2pt[iorder][iy][i][j] = new double[MaxPT]();
            for(int k=0;k<MaxPT;k++)
              dens2pt[iorder][iy][i][j][k] = 0.;
          }
      }
    }
  }
  // temporary
  double *** dens_tmp  = new double** [binRapidity];
  for(int iy=0;iy<binRapidity;iy++) {
    dens_tmp[iy] =  new double* [Maxx]();
    for(int i=0;i<Maxx;i++) {
        dens_tmp[iy][i] = new double[Maxy]();
        for (int j=0;j<Maxy;j++) dens_tmp[iy][i][j]=0;
    }
  }
  // temporary table which stores pt-unintegrated energy or entropy density
  double ****  dens_tmp_pt  = new double*** [binRapidity];
  for(int iy=0;iy<binRapidity;iy++) {
    dens_tmp_pt[iy] =  new double** [Maxx]();
    for(int i=0;i<Maxx;i++) {
        dens_tmp_pt[iy][i] = new double* [Maxy]();
        for (int j=0;j<Maxy;j++) {
          dens_tmp_pt[iy][i][j] = new double [MaxPT]();
          for(int k=0;k<MaxPT;k++)
            dens_tmp_pt[iy][i][j][k] = 0.;
        }
    }
  }
  // TA*TB profile (rotated using entropy):
  char file_TATB_Sd_4col[] = "data/TATB_fromSd_order_%d_4col.dat";
  char file_TATB_Sd_block[] = "data/TATB_fromSd_order_%d_block.dat";
  double **** TATB_Sd  = new double*** [number_of_orders];
  for(int iorder=0; iorder<number_of_orders; iorder++) // iorder starts from 0
  {
    TATB_Sd[iorder] = new double** [binRapidity];
    for(int iy=0;iy<binRapidity;iy++) {
      TATB_Sd[iorder][iy] =  new double* [Maxx]();
      for(int i=0;i<Maxx;i++) {
          TATB_Sd[iorder][iy][i] = new double[Maxy]();
          for (int j=0;j<Maxy;j++) TATB_Sd[iorder][iy][i][j]=0;
      }
    }
  }

  // TA*TB profile (rotated using energy):
  char file_TATB_Ed_4col[] = "data/TATB_fromEd_order_%d_4col.dat";
  char file_TATB_Ed_block[] = "data/TATB_fromEd_order_%d_block.dat";
  double **** TATB_Ed  = new double*** [number_of_orders];
  for(int iorder=0; iorder<number_of_orders; iorder++) // iorder starts from 0
  {
    TATB_Ed[iorder] = new double** [binRapidity];
    for(int iy=0;iy<binRapidity;iy++) {
      TATB_Ed[iorder][iy] =  new double* [Maxx]();
      for(int i=0;i<Maxx;i++) {
          TATB_Ed[iorder][iy][i] = new double[Maxy]();
          for (int j=0;j<Maxy;j++) TATB_Ed[iorder][iy][i][j]=0;
      }
    }
  }

  int output_TATB = paraRdr->getVal("output_TATB");

  long event = 1;

  // prepare to use ArraySaver class to auto-backup intermediate results
  long auto_backup_after_number_of_averaging = paraRdr->getVal("backup_number")-1;
  ArraySaver<double> dens1_saver("backup/dens1.dat", dens1, 4, Maxy, Maxx, binRapidity, number_of_orders);
  ArraySaver<double> dens2_saver("backup/dens2.dat", dens2, 4, Maxy, Maxx, binRapidity, number_of_orders);
  long backup_counter = auto_backup_after_number_of_averaging;
  ArraySaver<long> event_index_saver("backup/event_index.dat", &event, 1, 1);

  // event start.
  while (event<=nevent)
  {
    int tries = 0;
    int binary = 0;

    while (binary==0 || mc->CentralityCut()==0)
    {
      double b = sqrt((bmax*bmax - bmin*bmin)*drand48()+bmin*bmin);
      mc->generateNucleus(b,proj,targ);
      binary = mc->getBinaryCollision();
      //mc->dumpBinaryTable(); // for debugging
      if(binary==0 || mc->CentralityCut()==0) mc->deleteNucleus();
      tries++;
      if (tries>300)
      {
          cout << "===== MakeDensity::generate error =====" << endl
              << "No collisions detected after maximum number of tries." << endl
              << "What impact parameter are you using?" << endl;
          //exit(-1); // no more tries: something must be wrong.
          tries=0;
      }
    }
    Npart = mc->getNpart1()+mc->getNpart2();
    // need to rotate compute density before rotation.
    mc->getTA2();

    bool cutdSdypassFlag = true;
    for(int iorder=0; iorder<number_of_orders; iorder++)
    {
      int order = iorder + average_from_order; // "real" order used for averaging
      if(PTinte > 0)  // for integrated-pt
      {
        for(int iy=0;iy<binRapidity;iy++) {
          mc->setDensity(iy, -1);
          // cut total entropy density
          if(iy == 0 && cutdSdy == 1)
          {
             double totaldSdy = gettotaldSdy(iy);
             if(totaldSdy < cutdSdy_lowerBound || totaldSdy > cutdSdy_upperBound)
             {
                cutdSdypassFlag = false;
                break;
             }
          }
          // average entropy profile
          if (use_sd)
          {
              mc->rotateGrid(iy, order); // rotate grid according to gluon density <-> according to entropy density. Note that different rapidity slices are rotated separately, and this does not quite make sense.
              mc->getTA2();
              mc->setDensity(iy, -1); // now it's after rotation
              setSd(dens_tmp, iy); // includes factor multiplication
              // averaging --- entropy density:
              for(int i=0;i<Maxx;i++)
              for(int j=0;j<Maxy;j++)
              {
                  dens1[iorder][iy][i][j] = (dens1[iorder][iy][i][j]*(event-1) + dens_tmp[iy][i][j])/(double)(event); // event = number of succeeded events
              }
              // dumping TA*TB
              if (output_TATB) {
                  for(int i=0;i<Maxx;i++)
                  for(int j=0;j<Maxy;j++)
                  {
                      TATB_Sd[iorder][iy][i][j] = (TATB_Sd[iorder][iy][i][j]*(event-1) + mc->getTA1(i,j)*mc->getTA2(i,j))/(double)(event); // event = number of succeeded events
                  }
              }
          }
          // average energy profile
          if (use_ed)
          {
              setEd(dens_tmp, iy); // get energy density first
              // write energy density "back" to the "gluon density" matrix in MCnucl
              for(int i=0;i<Maxx;i++)
              for(int j=0;j<Maxy;j++) {
                  mc->setRho(iy,i,j,dens_tmp[iy][i][j]);
              }
              mc->rotateGrid(iy, order); // rotate grid according to energy density. Note that different rapidity slices are rotated separately, and this does not quite make sense.
              mc->getTA2();
              mc->setDensity(iy, -1); // now it's after rotation
              setEd(dens_tmp, iy); // includes factor multiplication
              // averaging --- entropy density:
              for(int i=0;i<Maxx;i++)
              for(int j=0;j<Maxy;j++)
              {
                  dens2[iorder][iy][i][j] = (dens2[iorder][iy][i][j]*(event-1) + dens_tmp[iy][i][j])/(double)(event); // event = number of succeeded events
              }
              // dumping TA*TB
              if (output_TATB) {
                  for(int i=0;i<Maxx;i++)
                  for(int j=0;j<Maxy;j++)
                  {
                      TATB_Ed[iorder][iy][i][j] = (TATB_Ed[iorder][iy][i][j]*(event-1) + mc->getTA1(i,j)*mc->getTA2(i,j))/(double)(event); // event = number of succeeded events
                  }
              }
          } // <->  if (use_ed)
        } // <-> for(int iy=0;iy<binRapidity;iy++)
      } //<-> if PTinte<0

      if(PTinte < 0)//pt is not integrated out
      {
        for(int iy=0;iy<binRapidity;iy++) 
        {
          mc->setDensity(iy,-1);  //generate pt-integrated table for grid rotation
          // cut total entropy density
          if(iy == 0 && cutdSdy == 1)
          {
             double totaldSdy = gettotaldSdy(iy);
             if(totaldSdy < cutdSdy_lowerBound || totaldSdy > cutdSdy_upperBound)
             {
                cutdSdypassFlag = false;
                break;
             }
          }
          //get dNdyd^2rdPt table
          for(int ipt=0;ipt<MaxPT;ipt++)
            mc->setDensity(iy,ipt);  

          // average entropy profile
          if (use_sd)
          {
            //Rotate according to the weighted center of pt-integrated dN/dyd^2r table
            mc->rotateGrid(iy, order); // rotate grid according to gluon density <-> according to entropy density. Note that different rapidity slices are rotated separately, and this does not quite make sense.
            mc->getTA2();

            for(int ipt = 0; ipt < MaxPT; ipt++)
              mc->setDensity(iy,ipt); // now it's after rotation
            setSd(dens_tmp_pt, iy); // includes factor multiplication
            // averaging --- entropy density:
            for(int i=0;i<Maxx;i++)
            for(int j=0;j<Maxy;j++)
            for(int ipt=0;ipt<MaxPT;ipt++)
            {
              dens1pt[iorder][iy][i][j][ipt] = (dens1pt[iorder][iy][i][j][ipt]*(event-1) + dens_tmp_pt[iy][i][j][ipt])/(double)(event); // event = number of succeeded events
            }
          }
            // average energy profile, not applicable
            if (use_ed)
            {
              //get weighted center using energy density as the weighting function
              setEd(dens_tmp, iy); // get energy density first
              // write energy density "back" to the "gluon density" matrix in MCnucl
              for(int i=0;i<Maxx;i++)
              for(int j=0;j<Maxy;j++) {
                mc->setRho(iy,i,j, dens_tmp[iy][i][j]);
              }
              mc->rotateGrid(iy, order); // rotate grid according to energy density. Note that different rapidity slices are rotated separately, and this does not quite make sense.
              mc->getTA2();

              for(int ipt = 0; ipt < MaxPT; ipt++)
                mc->setDensity(iy,ipt); // now it's after rotation
              setEd(dens_tmp_pt, iy); // includes factor multiplication
              // averaging --- entropy density:
              for(int i=0;i<Maxx;i++)
              for(int j=0;j<Maxy;j++)
              for(int ipt = 0; ipt < MaxPT; ipt++)
              {
                dens2pt[iorder][iy][i][j][ipt] = (dens2pt[iorder][iy][i][j][ipt]*(event-1) + dens_tmp_pt[iy][i][j][ipt])/(double)(event); // event = number of succeeded events
              }
            }
        } // <-> for(int iy=0;iy<binRapidity;iy++)
      }  //<-> if PTinte<0
      if(not cutdSdypassFlag)
        break;
    } // <-> for(int iorder=0; iorder<number_of_orders; iorder++)
    mc->deleteNucleus();

    // auto-backup
    if (backup_counter>0) backup_counter--; // autobackup turned on; keep track of counter
    else if (backup_counter==0) // autobackup turned on; time for output!
    {
      if (use_sd) dens1_saver.snapshot();
      if (use_ed) dens2_saver.snapshot();
      event_index_saver.snapshot();
      backup_counter = auto_backup_after_number_of_averaging;
    }

    if(cutdSdypassFlag)
    {
      cout << "processing event: " << event << endl;
      event++;
    }
  } // <-> while(event<events)

  // output average results
  char buffer[200];
  for (int iorder=0; iorder<number_of_orders; iorder++)
  {
    int order = iorder + average_from_order;
    for (int iy=0; iy<binRapidity; iy++)
    {
        // entropy
        if (use_sd) {
            if (use_4col)
            {
              sprintf(buffer, file1_4col, order);
              dumpDensity4Col(buffer, dens1[iorder], iy);
            }
            if (use_block)
            {
              sprintf(buffer, file1_block, order);
              dumpDensityBlock(buffer, dens1[iorder], iy);
            }
            if (use_5col)
            {
              sprintf(buffer, file1_5col, order);
              dumpDensity5Col(buffer, dens1pt[iorder], iy);
            }
        }
        // energy
        if (use_ed) {
            if (use_4col)
            {
              sprintf(buffer, file2_4col, order);
              dumpDensity4Col(buffer, dens2[iorder], iy);
            }
            if (use_block)
            {
              sprintf(buffer, file2_block, order);
              dumpDensityBlock(buffer, dens2[iorder], iy);
            }
           if (use_5col)
            {
              sprintf(buffer, file2_5col, order);
              dumpDensity5Col(buffer, dens2pt[iorder], iy);
            }            
        }
        // output TA*TB
        if (output_TATB) {
            // TA*TB from entropy
            if (use_sd) {
                if (use_4col)
                {
                  sprintf(buffer, file_TATB_Sd_4col, order);
                  dumpDensity4Col(buffer, TATB_Sd[iorder], iy);
                }
                if (use_block)
                {
                  sprintf(buffer, file_TATB_Sd_block, order);
                  dumpDensityBlock(buffer, TATB_Sd[iorder], iy);
                }
            }
            // TA*TB from energy
            if (use_ed) {
                if (use_4col)
                {
                  sprintf(buffer, file_TATB_Ed_4col, order);
                  dumpDensity4Col(buffer, TATB_Ed[iorder], iy);
                }
                if (use_block)
                {
                  sprintf(buffer, file_TATB_Ed_block, order);
                  dumpDensityBlock(buffer, TATB_Ed[iorder], iy);
                }
            }
        }
    }
  }

  // clean up
  for(int iorder=0; iorder<number_of_orders; iorder++) {
    for(int iy=0;iy<binRapidity;iy++) {
      for(int i=0;i<Maxx;i++) delete [] dens1[iorder][iy][i];
      delete [] dens1[iorder][iy];
    }
    delete [] dens1[iorder];
  }
  delete [] dens1;

  for(int iorder=0; iorder<number_of_orders; iorder++) {
    for(int iy=0;iy<binRapidity;iy++) {
      for(int i=0;i<Maxx;i++) delete [] dens2[iorder][iy][i];
      delete [] dens2[iorder][iy];
    }
    delete [] dens2[iorder];
  }
  delete [] dens2;


  for(int iy=0;iy<binRapidity;iy++) {
    for(int i=0;i<Maxx;i++) {
      for(int j=0;j<Maxy;j++) delete [] dens_tmp_pt[iy][i][j];
        delete [] dens_tmp_pt[iy][i];
    }
    delete [] dens_tmp_pt[iy];
  }
  delete [] dens_tmp_pt;

  for(int iorder=0; iorder<number_of_orders; iorder++) {
    for(int iy=0;iy<binRapidity;iy++) {
      for(int i=0;i<Maxx;i++) {
        for(int j=0;j<Maxy;j++) {
          delete [] dens1pt[iorder][iy][i][j];
          delete [] dens2pt[iorder][iy][i][j];
        }
        delete [] dens1pt[iorder][iy][i];
        delete [] dens2pt[iorder][iy][i];
      }
      delete [] dens1pt[iorder][iy];
      delete [] dens2pt[iorder][iy];
    }
    delete [] dens1pt[iorder];
    delete [] dens2pt[iorder];
  }
  delete [] dens1pt;
  delete [] dens2pt;

  for(int iy=0;iy<binRapidity;iy++) {
    for(int i=0;i<Maxx;i++) delete [] dens_tmp[iy][i];
    delete [] dens_tmp[iy];
  }
  delete [] dens_tmp;

  for(int iorder=0; iorder<number_of_orders; iorder++) {
    for(int iy=0;iy<binRapidity;iy++) {
      for(int i=0;i<Maxx;i++) delete [] TATB_Sd[iorder][iy][i];
      delete [] TATB_Sd[iorder][iy];
    }
    delete [] TATB_Sd[iorder];
  }
  delete [] TATB_Sd;

  for(int iorder=0; iorder<number_of_orders; iorder++) {
    for(int iy=0;iy<binRapidity;iy++) {
      for(int i=0;i<Maxx;i++) delete [] TATB_Ed[iorder][iy][i];
      delete [] TATB_Ed[iorder][iy];
    }
    delete [] TATB_Ed[iorder];
  }
  delete [] TATB_Ed;

}
//----------------------------------------------------------------------


//**********************************************************************
void MakeDensity::generateEccTable(int nevent)
// bmin, bmax: possible impact parameter range to be sampled
// Npmin, Npmax: cut on N_part
{
  int from_order=paraRdr->getVal("ecc_from_order"), to_order=paraRdr->getVal("ecc_to_order");

  // entropy profile:
  char file1[] = "data/sn_ecc_eccp_%d.dat";
  double *** dens1  = new double** [binRapidity];
  for(int iy=0;iy<binRapidity;iy++) {
    dens1[iy] =  new double* [Maxx]();
    for(int i=0;i<Maxx;i++) {
        dens1[iy][i] = new double[Maxy]();
        for (int j=0;j<Maxy;j++) dens1[iy][i][j]=0;
    }
  }

  // energy profile:
  char file2[] = "data/en_ecc_eccp_%d.dat";
  double *** dens2  = new double** [binRapidity];
  for(int iy=0;iy<binRapidity;iy++) {
    dens2[iy] =  new double* [Maxx]();
    for(int i=0;i<Maxx;i++) {
        dens2[iy][i] = new double[Maxy];
        for (int j=0;j<Maxy;j++) dens2[iy][i][j]=0;
    }
  }

  double b;

  // event start.
  int event=1;
  while (event<=nevent)
  {
    int tries = 0;
    int binary = 0;
    while (binary==0 || mc->CentralityCut()==0)
    {
      b = sqrt((bmax*bmax - bmin*bmin)*drand48()+bmin*bmin);
      mc->generateNucleus(b,proj,targ);
      binary = mc->getBinaryCollision();
      if(binary==0 || mc->CentralityCut()==0) mc->deleteNucleus();
      tries++;
      if (tries>300)
      {
          cout << "===== MakeDensity::generate error =====" << endl
              << "No collisions detected after maximum number of tries." << endl
              << "What impact parameter are you using?" << endl;
          //exit(-1); // no more tries: something must be wrong.
          tries=0;
      }
    }
    Npart = mc->getNpart1()+mc->getNpart2();
    Nbin = mc->getNcoll();

    // compute eccentricity.
    mc->getTA2();
    
    bool cutdSdypassFlag = true;
    for(int iy=0;iy<binRapidity;iy++) {
      mc->setDensity(iy, -1);
      // cut total entropy
      if(iy == 0 && cutdSdy == 1)
      {
         double totaldSdy = gettotaldSdy(iy);
         if(totaldSdy < cutdSdy_lowerBound || totaldSdy > cutdSdy_upperBound)
         {
            cutdSdypassFlag = false;
            break;
         }
      }
      // entropy first
      if (paraRdr->getVal("use_sd"))
      {
        setSd(dens1, iy);
        dumpEccentricities(file1, dens1, iy, from_order, to_order, Npart, Nbin, b);
      }
      // energy second
      if (paraRdr->getVal("use_ed"))
      {
        setEd(dens2, iy);
        dumpEccentricities(file2, dens2, iy, from_order, to_order, Npart, Nbin, b);
      }
    }

    mc->deleteNucleus();
    if(cutdSdypassFlag)
    {
      cout << "processing event: " << event << endl;
      event++;
    }
  } // <-> while (event<=nevent)


  for(int iy=0;iy<binRapidity;iy++) {
    for(int j=0;j<Maxy;j++) delete [] dens1[iy][j];
    delete [] dens1[iy];
  }
  delete [] dens1;

  for(int iy=0;iy<binRapidity;iy++) {
    for(int j=0;j<Maxy;j++) delete [] dens2[iy][j];
    delete [] dens2[iy];
  }
  delete [] dens2;

}
//----------------------------------------------------------------------

//********************************************************************
void MakeDensity::dumpEccentricities(char* base_filename, double*** density, const int iy, int from_order, int to_order, double Npart_current, double Nbin_current, double b)
// calculate and output eccentricities
{
    bool deformedFlag = false;
    if(paraRdr->getVal("proj_deformed") == 1 or paraRdr->getVal("targ_deformed") == 1) 
      deformedFlag = true;
    std::ofstream of;
    char buffer[200];
    double x,y,r,theta;
    int order;
    double  xc, yc, total,
            *mom_real = new double[to_order+1], *mom_imag = new double[to_order+1], *norm = new double[to_order+1], // norm is <r^2>
            *momp_real = new double[to_order+1], *momp_imag = new double[to_order+1], *normp = new double[to_order+1]; // normp is <r^n>

    // center of mass:
    xc = 0.0; yc = 0.0; total = 0.0;
    for(int i=0;i<Maxx;i++)
    for(int j=0;j<Maxy;j++)
    {
        x = Xmin + i*dx; y = Ymin + j*dy;
        xc += x*density[iy][i][j];
        yc += y*density[iy][i][j];
        total += density[iy][i][j];
    }
    xc /= total; yc /= total;

/*
    // for overlapping area:
    // method 1, pi*sqrt(<x^2> <y^2>)
    double x2 = 0.0, y2 = 0.0, overlap_area1 = 0.0; // <x^2>, <y^2>, and different definitions for overlapping area
    for(int i=0;i<Maxx;i++)
    for(int j=0;j<Maxy;j++)
    {
        x = Xmin + i*dx - xc; y = Ymin + j*dy - yc;
        x2 += x*x*density[iy][i][j];
        y2 += y*y*density[iy][i][j];
    }
    x2 /= total; y2 /= total;
    overlap_area1 = M_PI*sqrt(x2*y2);
    // method 0, use the area inside which the integral is 0.393469
    vector<double> density_ordered_vec(Maxx*Maxy);
    for(int i=0;i<Maxx;i++)
    for(int j=0;j<Maxy;j++)
    {
        density_ordered_vec[i*Maxy+j] = density[iy][i][j] / total;
    }
    // now sort it and calculate the area inside which the integral is 0.39
    sort(density_ordered_vec.begin(), density_ordered_vec.end());
    double overlap_sum = 0.0, overlap_area2 = 0.0;
    for (long i=0; i<Maxx*Maxy; i++)
    {
        overlap_sum += density_ordered_vec[Maxx*Maxy-1-i];
        if (overlap_sum > 0.393469)
        {
            //cout << i << "   " << density_ordered_vec[i] << endl;
            overlap_area2 = i*dx*dy;
            break;
        }
    }
    //cout << overlap_sum << "   " << overlap_area1 << "   " << overlap_area2 << endl;
*/
    // for eccentricity:
    for (order=from_order; order<=to_order; order++)
    {
        // initialize:
        mom_real[order]=0.0; mom_imag[order]=0.0; norm[order]=0.0;
        momp_real[order]=0.0; momp_imag[order]=0.0; normp[order]=0.0;

        // calculate numerator and denominators:
        for(int i=0;i<Maxx;i++)
        for(int j=0;j<Maxy;j++)
        {
            x = Xmin + i*dx - xc; y = Ymin + j*dy - yc; // shift to center
            r = sqrt(x*x + y*y); theta = atan2(y,x);
            mom_real[order] += r*r*cos(order*theta)*density[iy][i][j];
            mom_imag[order] += r*r*sin(order*theta)*density[iy][i][j];
            norm[order] += r*r*density[iy][i][j];
            if(order == 1)
            {
               momp_real[order] += pow(r,3)*cos(order*theta)*density[iy][i][j];
               momp_imag[order] += pow(r,3)*sin(order*theta)*density[iy][i][j];
               normp[order] += pow(r,3)*density[iy][i][j];
            }
            else
            {
               momp_real[order] += pow(r,order)*cos(order*theta)*density[iy][i][j];
               momp_imag[order] += pow(r,order)*sin(order*theta)*density[iy][i][j];
               normp[order] += pow(r,order)*density[iy][i][j];
            }
        }

        // take ratio; note that the minus sign is just a convention
        mom_real[order] = -mom_real[order]/norm[order];
        mom_imag[order] = -mom_imag[order]/norm[order];
        momp_real[order] = -momp_real[order]/normp[order];
        momp_imag[order] = -momp_imag[order]/normp[order];
        // and output:
        sprintf(buffer, base_filename, order);
        of.open(buffer, std::ios_base::app);
        
        if(deformedFlag)
        {
            of  << setprecision(8) << setw(16) <<  mom_real[order]
                << setprecision(8) << setw(16) <<  mom_imag[order]
                << setprecision(8) << setw(16) <<  momp_real[order]
                << setprecision(8) << setw(16) <<  momp_imag[order]
                << setprecision(5) << setw(10)  <<  Npart_current
                << setprecision(5) << setw(10) <<  Nbin_current
                << setprecision(8) << setw(16) <<  total*dx*dy // integrated profile measure
                << setprecision(8) << setw(16) <<  b
//                << setprecision(8) << setw(16) <<  overlap_area1
//                << setprecision(8) << setw(16) <<  overlap_area2
                << setprecision(8) << setw(16) <<  mc->lastCx1
                << setprecision(8) << setw(16) <<  mc->lastPh1
                << setprecision(8) << setw(16) <<  mc->lastCx2
	          << setprecision(8) << setw(16) <<  mc->lastPh2
//                << setprecision(8) << setw(16) <<  xc
//                << setprecision(8) << setw(16) <<  yc
//                << setprecision(3) << setw(12) <<  rapMin+(rapMax-rapMin)/binRapidity*iy
                << endl;
        }
        else
        {
            of  << setprecision(8) << setw(16) <<  mom_real[order]
                << setprecision(8) << setw(16) <<  mom_imag[order]
                << setprecision(8) << setw(16) <<  momp_real[order]
                << setprecision(8) << setw(16) <<  momp_imag[order]
                << setprecision(5) << setw(10) <<  Npart_current
                << setprecision(5) << setw(10) <<  Nbin_current
                << setprecision(8) << setw(16) <<  total*dx*dy // integrated profile measure
                << setprecision(8) << setw(16) <<  b
                << endl;
        }
        of.close();
    }

    delete[] mom_real, mom_imag, norm, momp_real, momp_imag, normp;
}


void MakeDensity::setSd(double*** data, const int iy)
// Copy density from mc object to data using mc->getRho function, but with multiplicity factor (gluon density to entropy)
{
    for(int i=0;i<Maxx;i++)
    for(int j=0;j<Maxy;j++)
    {
        data[iy][i][j] = mc->getRho(iy,i,j)*finalFactor;
    }
}

void MakeDensity::setEd(double*** data, const int iy)
// Copy density from mc object to data using mc->getRho function, with multiplicity factor, and converts to energy
{
    for(int i=0;i<Maxx;i++)
    for(int j=0;j<Maxy;j++)
    {
        data[iy][i][j] = eos.edFromSd(mc->getRho(iy,i,j)*finalFactor);
    }
}

void MakeDensity::setSd(double**** data, const int iy)
  // Copy density from mc object to data using mc->getRho function, but with multiplicity factor (gluon density to entropy)
{
    for(int i=0;i<Maxx;i++)
    for(int j=0;j<Maxy;j++)
    for(int ipt=0;ipt<MaxPT;ipt++)
    {
        data[iy][i][j][ipt] = mc->getRho(iy,i,j,ipt)*finalFactor;
    }
}

void MakeDensity::setEd(double**** data, const int iy)
// Copy density from mc object to data using mc->getRho function, with multiplicity factor, and converts to energy
{
    for(int i=0;i<Maxx;i++)
    for(int j=0;j<Maxy;j++)
    for(int ipt=0;ipt<MaxPT;ipt++)
    {
        data[iy][i][j][ipt] = eos.edFromSd(mc->getRho(iy,i,j,ipt)*finalFactor);
    }
}

double MakeDensity::gettotaldSdy(const int iy)
// calculate total entropy density at rapidity iy using mc->getRho function with multiplicity factor (gluon density to entropy)
{
    double sum = 0.0;
    for(int i=0;i<Maxx;i++)
    for(int j=0;j<Maxy;j++)
    {
        sum += mc->getRho(iy,i,j)*dx*dy;
    }
    return(sum);
}

//--------------------------------------------------------------------


//********************************************************************
void MakeDensity::dumpDensity4Col(char filename[], double *** data, const int iy)
{
  ofstream of;
  of.open(filename, std::ios_base::app);
  of  << "# <npart>= " << Npart
      << " xmax= " << Maxx << " ymax= " << Maxy
      << endl;

    for(int i=0;i<Maxx;i++)
    for(int j=0;j<Maxy;j++)
    {
      double x = Xmin + i*dx;
      double y = Ymin + j*dy;
      double rap = rapMin+(rapMax-rapMin)/binRapidity*iy;
      of <<  setprecision(3) << setw(10) <<  rap
          << setprecision(3) << setw(10) <<  x
          << setprecision(3) << setw(10) <<  y
          << setprecision(12) << setw(22) <<  data[iy][i][j]
	  << endl;
    }
  of.close();
}

void MakeDensity::dumpDensityBlock(char filename[], double *** data, const int iy)
{
  ofstream of;
  of.open(filename, std::ios_base::app);
  for(int i=0;i<Maxx;i++) {
    for(int j=0;j<Maxy;j++) {
      double x = Xmin + i*dx;
      double y = Ymin + j*dy;
      of << setprecision(12) << setw(22) << data[iy][i][j];
    }
    of << endl;
  }
  of.close();
}

void MakeDensity::dumpDensity5Col(char filename[], double **** data, const int iy)
{
  ofstream of;
  of.open(filename, std::ios_base::out);
  of  << "% <npart>= " << Npart
      << " xmax= " << Maxx << " ymax= " << Maxy
      << " Ptmax= " << MaxPT
      << endl;

    for(int i=0;i<Maxx;i++)
    for(int j=0;j<Maxy;j++)
    {
      for(int ipt=0;ipt<MaxPT;ipt++)
      {
        // double x = Xmin + i*dx;
        // double y = Ymin + j*dy;
        // double ptstep = PTmin + ipt*dpt;
        // double rap = rapMin+(rapMax-rapMin)/binRapidity*iy;
        // of  << setprecision(3) << setw(14) <<  rap
        //     << setprecision(3) << setw(14) <<  x
        //     << setprecision(3) << setw(14) <<  y
        //     << setprecision(3) << setw(14) <<  ptstep;
        of  << scientific << setprecision(12) << setw(22) <<  data[iy][i][j][ipt];
      }
      of << endl;
    }

  of.close();
}

void MakeDensity::dumpDesityptCol(char filename[], double **** data, const int iy)
{
// a concise way of storing pt-differential table 
  ofstream of;
  of.open(filename, std::ios_base::out);
  of  << "% <npart>= " << Npart
      << " xmax= " << Maxx << " ymax= " << Maxy
      << " Ptmax= " << MaxPT
      << endl;

  for(int i=0;i<Maxx;i++)
    for(int j=0;j<Maxy;j++)
    {
      for(int ipt=0;ipt<MaxPT;ipt++)
        of  << scientific << setprecision(12) << setw(22) <<  data[iy][i][j][ipt];
      of << endl;
    }

  of.close();
}

//--------------------------------------------------------------------
