#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <string>
#include <cstdlib>

#include "MathBasics.h"
#include "OverLap.h"
#include "ParameterReader.h"

using namespace std;

//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
// --- initializes parameters for WS/Hulthen density distribution of a nucleus
//     function getRandomWS(x,y,z) then returns sampled nucleon coordinate
OverLap::OverLap(ParameterReader* paraRdr_in, int a, double signnin, int deformed_in)
{
  paraRdr = paraRdr_in;

  // save atomic # of nucleus for later; can be recalled using: getAtomic()
  deformed = deformed_in;
  A = (double)a;
  atomic = a;

  // if nucleon, nothing more needed
  if (a==1) return;

  // generic (crude) default parameterization of radius and surface thickness
  rad=1.12*pow(A,0.333333)-0.86/pow(A,0.333333);
  dr = 0.54;
  density0 = 0.17;

  //taken from C.W.De Jager et al. Atom.Data.Nucl.Data Tabl.36, 495 (1987).
  //if(a==197){
  //  rad = 6.38;
  //  dr = 0.535;
  //  density0 = 0.1695;
  //}else if(a==63){
  //  rad = 4.20641;
  //  dr = 0.5877;
  //  density0 = 0.1686;
  //}

  // to reproduce ordinary Woods-Saxon for finite size nucleon (T. Hirano)
  if(a == 197){
    //Reparametrization (Au)
    rad=6.42;
    dr=0.45;
    density0 = 0.1695;
  }else if(a==63){
    //Reparametrization (Cu)
    rad=4.28;
    dr=0.5;
    density0 = 0.1686;
  }else if(a==238){
    //Reparametrization (U)
    rad=6.86;
    dr=0.44;
    density0 = 0.166;

    //Taken from P.Filip et al.,PRC80,054903(2009)
    //      rad=6.81;
    //      dr=0.54;
    //      density0 = 0.17;
  }else if(a==208){
    //Reparametrization (Pb)
    rad = 6.67;
    dr = 0.44;
    density0 = 0.161;
  }

  // default deformation parameters, 0: no deform, 1: deform (05/03/2010 by TH)
  beta2 = 0.0;
  beta4 = 0.0;
  ctr = 1.0;
  phir = 0.0;
  rmaxCut = rad + 2.5;
  rwMax = 1.0/(1.0+exp(-rad/dr));

  sig = signnin/10.0;
  zini = -10.0;
  zfin = 10.0;
  Gauss38(zini,zfin,z,zw);

  if(deformed){
    // Parameters taken from Moller et al., Atomic Data and Nuclear Data
    // Tables, 59, 185-381, (1995).
    if(atomic == 197){
      //(Au)
      beta2 = -0.13;
      beta4 = -0.03;
    }else if(atomic == 63){
      //(Cu)
      beta2 = 0.162;
      beta4 = 0.006;
    }else if(atomic == 238){
      //(U)
      beta2 = 0.28;
      beta4 = 0.093;
    }else{
      cerr << "Mass number not available for reparametrization" << endl;
    }
  }
  else
  {
    beta2 = 0.0; beta4 = 0.0;
  }

  if (atomic == 3) // read in triton position
  {
    readin_triton_position();
  }

  flag_NN_correlation = paraRdr->getVal("include_NN_correlation");
  if(flag_NN_correlation == 1 && (atomic == 197 || atomic == 208))
  {
     readin_nucleon_positions();
  }
  
}

OverLap::~OverLap()
{
   if(flag_NN_correlation == 1 && (atomic == 197 || atomic == 208))
   {
      for(int iconf = 0; iconf < n_configuration; iconf++)
      {
         for(int ia = 0; ia < atomic; ia++)
            delete [] nucleon_pos_array[iconf][ia];
         delete [] nucleon_pos_array[iconf];
      }
      delete [] nucleon_pos_array;
   }
}

void OverLap::GetDeuteronPosition(double& x1,double& y1,double& z1,double& x2,double& y2,double& z2)
{
   //get proton/neutron separation d (fm)
   double d;
   
   d=sample_deuteron.rand(); //now sample random rotation of deuteron
   x1 = (d/2.0);
   y1 = 0;
   z1 = 0;

   Point3D p3d(x1,y1,z1);
   p3d.rotate(ctr, phir);
   x1 = p3d.x; y1 = p3d.y; z1 = p3d.z;

   x2 = -x1;
   y2 = -y1;
   z2 = -z1;
}

void OverLap::readin_triton_position()
{
   ifstream triton_position("tables/triton_positions.dat");
   double x1, y1, z1, x2, y2, z2, x3, y3, z3;
   triton_position >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> x3 >> y3 >> z3;
   while(!triton_position.eof())
   {
       vector<double> temp;
       temp.push_back(x1);
       temp.push_back(y1);
       temp.push_back(z1);
       temp.push_back(x2);
       temp.push_back(y2);
       temp.push_back(z2);
       temp.push_back(x3);
       temp.push_back(y3);
       temp.push_back(z3);
       triton_pos.push_back(temp);
       triton_position >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> x3 >> y3 >> z3;
   }
}

void OverLap::readin_nucleon_positions()
{
   cout << "read in nucleon positions with nucleon-nucleon correlation..." << flush;
   ostringstream filename;
   if(atomic == 197)
   {
      filename << "tables/au197-sw-full_3Bchains-conf1820.dat";
      n_configuration = 1820;
   }
   else if(atomic == 208)
   {
      //int temp = rand() % 10 + 1;
      int temp = 1;
      filename << "tables/pb208-" << temp << ".dat";
      n_configuration = 10000;
   }
   else
      return;

   nucleon_pos_array = new double** [n_configuration];
   for(int iconf = 0; iconf < n_configuration; iconf++)
   {
      nucleon_pos_array[iconf] = new double* [atomic];
      for(int ia = 0; ia < atomic; ia++)
         nucleon_pos_array[iconf][ia] = new double [3];
   }

   ifstream input(filename.str().c_str());
   for(int iconf = 0; iconf < n_configuration; iconf++)
   {
      for(int ia = 0; ia < atomic; ia++)
      {
         double x_local, y_local, z_local;
         int isospin;
         int dummy;
         if(atomic == 208)
            input >> x_local >> y_local >> z_local >> isospin;
         else
            input >> x_local >> y_local >> z_local >> isospin >> dummy;
         nucleon_pos_array[iconf][ia][0] = x_local;
         nucleon_pos_array[iconf][ia][1] = y_local;
         nucleon_pos_array[iconf][ia][2] = z_local;
      }
   }
   input.close();
   cout << " done." << endl;
}

void OverLap::GetTritonPosition(double& x1,double& y1,double& z1,double &x2,double& y2,double& z2, double &x3, double &y3, double &z3)
{
   int num_configuration = triton_pos.size();
   int rand_num = rand() % num_configuration;

   x1 = triton_pos[rand_num][0];
   y1 = triton_pos[rand_num][1];
   z1 = triton_pos[rand_num][2];
   x2 = triton_pos[rand_num][3];
   y2 = triton_pos[rand_num][4];
   z2 = triton_pos[rand_num][5];
   x3 = triton_pos[rand_num][6];
   y3 = triton_pos[rand_num][7];
   z3 = triton_pos[rand_num][8];

   Point3D p3d1(x1,y1,z1);
   p3d1.rotate(ctr, phir);
   x1 = p3d1.x; y1 = p3d1.y; z1 = p3d1.z;

   Point3D p3d2(x2,y2,z2);
   p3d2.rotate(ctr, phir);
   x2 = p3d2.x; y2 = p3d2.y; z2 = p3d2.z;

   Point3D p3d3(x3,y3,z3);
   p3d3.rotate(ctr, phir);
   x3 = p3d3.x; y3 = p3d3.y; z3 = p3d3.z;
}

// *** this function applies to deformed nuclei ***
void OverLap::getDeformRandomWS(double& x, double& y, double& z)
{
  double rad1 = rad;
  double rwMax1 = rwMax;
  double r = 0.0;
  double cx = 1.0; // cx is cos(theta)

  if (deformed)
  {
    do {
      //Uniform distribution in a sphere with r = rmaxCut
      r = rmaxCut*pow(drand48(),1.0/3.0);
      cx = 1.0-2.0*drand48();

      //Deformation
      //Main axis in z-axis
      double y20 = SphericalHarmonics(2,cx);
      double y40 = SphericalHarmonics(4,cx);
      rad1 =rad*(1.0 + beta2 * y20 + beta4 * y40);
      rwMax1 =1.0/(1.0+exp(-rad1/dr));

    } while(drand48()*rwMax1 > 1.0/(1.0+exp((r-rad1)/dr)));

    double sx = sqrt(1.0-cx*cx); // sx is sin(theta)
    double phi = 2*M_PI*drand48();

    Point3D p3d(r*sx*cos(phi), r*sx*sin(phi), r*cx);
    p3d.rotate(ctr, phir);
    x = p3d.x; y = p3d.y; z = p3d.z;
  }
  else
  {
    do {
      //Uniform distribution in a sphere with r = rmaxCut
      r = rmaxCut*pow(drand48(),1.0/3.0);
    } while(drand48()*rwMax > 1.0/(1.0+exp((r-rad)/dr)));
    cx = 1.0-2.0*drand48();
    double sx = sqrt(1.0-cx*cx); // sx is sin(theta)
    double phi = 2*M_PI*drand48();
    x = r*sx*cos(phi);
    y = r*sx*sin(phi);
    z = r*cx;
  }
}

void OverLap::get_nucleon_position_with_NN_correlation(double **nucleon_ptr)
{
   int i_configuration = rand() % n_configuration;  // pick the configuration
   double** temp_nucleus = new double* [atomic];
   for(int ia = 0; ia < atomic; ia++)
      temp_nucleus[ia] = new double [3];

   // first recenter the nucleus
   double xcm = 0.0, ycm = 0.0, zcm = 0.0;
   for(int ia = 0; ia < atomic; ia++)
   {
      xcm += nucleon_pos_array[i_configuration][ia][0];
      ycm += nucleon_pos_array[i_configuration][ia][1];
      zcm += nucleon_pos_array[i_configuration][ia][2];
   }
   for(int ia = 0; ia < atomic; ia++)
   {
      temp_nucleus[ia][0] = nucleon_pos_array[i_configuration][ia][0] - xcm/atomic;
      temp_nucleus[ia][1] = nucleon_pos_array[i_configuration][ia][1] - ycm/atomic;
      temp_nucleus[ia][2] = nucleon_pos_array[i_configuration][ia][2] - zcm/atomic;
   }

   // then rotate a random angle
   double cos_theta = 1.0-2.0*drand48();
   double phi = 2*M_PI*drand48();
   for(int ia = 0; ia < atomic; ia++)
   {
      Point3D p3d(temp_nucleus[ia][0], temp_nucleus[ia][1], temp_nucleus[ia][2]);
      p3d.rotate(cos_theta, phi);
      temp_nucleus[ia][0] = p3d.x; 
      temp_nucleus[ia][1] = p3d.y; 
      temp_nucleus[ia][2] = p3d.z;
   }

   // return back
   for(int ia = 0; ia < atomic; ia++)
      for(int ii = 0; ii < 3; ii++)
         nucleon_ptr[ia][ii] = temp_nucleus[ia][ii];

   // clean up
   for(int ia = 0; ia < atomic; ia++)
      delete [] temp_nucleus[ia];
   delete [] temp_nucleus;
}


double OverLap::SphericalHarmonics(int l, double ct)
{
  //Currently assuming m=0 and available for Y_{20} and Y_{40}
  // "ct" is cos(theta)

  double ylm = 0.0;

  if(l == 2){

    ylm = 3.0*ct*ct-1.0;
    ylm *= 0.31539156525252005; //pow(5.0/16.0/M_PI,0.5);

  }else if (l == 4){

    ylm = 35.0*ct*ct*ct*ct;
    ylm -= 30.0*ct*ct;
    ylm += 3.0;
    ylm *= 0.10578554691520431; //3.0/16.0/pow(M_PI,0.5);

  }else{
    cerr << "Not available in Overlap::SphericalHarmonics" << endl;
  }

  return ylm;

}


//**********************************************************************
void OverLap::Gauss38(double xini,double xfin,double* xn,double* wn)
{
    double  x[38], w[38];

    x[37]=9.980499305357e-1;
    x[36]=9.897394542664e-1;
    x[35]=9.748463285902e-1;
    x[34]=9.534663309335e-1;
    x[33]=9.257413320486e-1;
    x[32]=8.918557390046e-1;
    x[31]=8.520350219324e-1;
    x[30]=8.065441676053e-1;
    x[29]=7.556859037540e-1;
    x[28]=6.997986803792e-1;
    x[27]=6.392544158297e-1;
    x[26]=5.744560210478e-1;
    x[25]=5.058347179279e-1;
    x[24]=4.338471694324e-1;
    x[23]=3.589724404794e-1;
    x[22]=2.817088097902e-1;
    x[21]=2.025704538921e-1;
    x[20]=1.220840253379e-1;
    x[19]=4.078514790458e-2;

//    .....   WEIGHT       ...........
    w[37]=5.002880749632e-3;
    w[36]=1.161344471647e-2;
    w[35]=1.815657770961e-2;
    w[34]=2.457973973823e-2;
    w[33]=3.083950054518e-2;
    w[32]=3.689408159400e-2;
    w[31]=4.270315850467e-2;
    w[30]=4.822806186076e-2;
    w[29]=5.343201991033e-2;
    w[28]=5.828039914700e-2;
    w[27]=6.274093339213e-2;
    w[26]=6.678393797914e-2;
    w[25]=7.038250706690e-2;
    w[24]=7.351269258474e-2;
    w[23]=7.615366354845e-2;
    w[22]=7.828784465821e-2;
    w[21]=7.990103324353e-2;
    w[20]=8.098249377060e-2;
    w[19]=8.152502928039e-2;

    for(int i=0;i<19;i++) {
      x[i] = -x[37-i];
      w[i] =  w[37-i];
    }
    for(int i=0;i<38;i++) {
      xn[i] =(xfin-xini)*x[i]/2.0+(xini+xfin)/2.0;
      wn[i] =(xfin-xini)*w[i]/2.0;
    }

}



//#define MAIN

#ifdef MAIN
#include <iomanip>
int main() {

  double r,x,y,z;
  int hist[100];
  double sig=42.0; // 200GeV.
  OverLap* overlap = new OverLap(2,0.0,sig,0);

  int count=0;
  for(int i=0; i<100; i++)  hist[i]=0;
  double Rrms=0.;

  for(int i=0;i<100000;i++) {
    overlap->getRandomWS(x,y,z);
    r = sqrt(x*x+y*y+z*z);
    hist[(int)(r/.1)]++;
    count++;
  }

  for(int i=0; i<100; i++) {
    Rrms += pow((i+.5)*.1,2.) * (float)hist[i]/(float)count;
    cout << (i+.5)*.1 << " " << (float)hist[i]/(float)count << endl;
  }
  cout << "# R_rms = " << sqrt(Rrms) << endl;

  return 0;
}
#endif
