#include <cmath>
#include <iomanip>
#include <vector>
#include <sstream>
#include "Freestreaming.h"
#include "gauss_quadrature.h"
#include "mistools.h"
#include "arsenal.h"

using namespace std;

FreeStrm::FreeStrm(double xmax, double ymax, double dx0,double dy0,
		    int ny, double rapmin, double rapmax, double taui, double tauf)
{
    Xmax=xmax;
    Ymax=ymax;
    Xmin=-xmax;
    Ymin=-ymax;

    dx=dx0;
    dy=dy0;
    Maxx=(int)((Xmax-Xmin)/dx+0.1)+1;
    Maxy=(int)((Ymax-Ymin)/dy+0.1)+1;
    nRap=ny;
    rapMin=rapmin;
    rapMax=rapmax;

    Taui=taui;
    Tauf=tauf;
    //for generating test profile
    PTmin = 0.1, dpt=0.1, PTmax=12.;
    MaxPT = (int)((PTmax-PTmin)/dpt + 0.1)+1;

    Xcm = 0.;
    Ycm = 0.;
    shiftedTable = 0;
    unshiftedTable = 0;
//  cout << "Free-streaming procedure Initialized!" << endl;
}


FreeStrm::~FreeStrm()
{  
  //clean up readin table
  if(unshiftedTable!=0)
  {
    for(int iy=0;iy<nRap;iy++)  {
        {
          for(int i=0;i<Maxx;i++) 
            delete [] unshiftedTable[iy][i];
        }
      delete [] unshiftedTable[iy];
      }
    delete [] unshiftedTable;      
  }    

  //Clean intermediate table
  if(shiftedTable!=0)
  {
    for(int iy=0;iy<nRap;iy++)  {
        {
          for(int i=0;i<Maxx;i++) 
            delete [] shiftedTable[iy][i];
        }
      delete [] shiftedTable[iy];
      }
    delete [] shiftedTable;      
  }  
}

void FreeStrm::CopyTable(double ***source)
{
  //store readin table
  unshiftedTable  = new double** [nRap];
  for(int iy=0;iy<nRap;iy++) {
  unshiftedTable[iy] =  new double* [Maxx];
  for(int i=0;i<Maxx;i++) {
      unshiftedTable[iy][i] = new double [Maxy];
      for(int j=0;j<Maxy;j++) {
          unshiftedTable[iy][i][j]=source[iy][i][j];
      }
    }
  }  
}

double FreeStrm::GetDensity(int iRap, int i, int j, double phip)
{
/* Get the shifted profile at (i,j) in x-y plane for a specific angle phip and travling 
   time DTau
*/
  double x,y;     //unshifted coordinate
  double shfedx, shfedy;    //shifed coordinate
  double DTau, Phip;
  double it, jt;             //indices corresponding to shifted coordinates

  Phip=phip;

  x=Xmin+i*dx;
  y=Ymin+j*dy;
  DTau=Tauf-Taui;

  shfedx=x-DTau*cos(Phip);
  shfedy=y-DTau*sin(Phip);


  if( (shfedx>Xmax||shfedx<Xmin) ||(shfedy>Ymax||shfedy<Ymin) )  //coordinates locate outside of the grid
      return 0;

  else
  {   
    it=(shfedx-Xmin)/dx;        
    jt=(shfedy-Ymin)/dy;

    // boundary safty control
    if (jt<=0) jt += 1e-10;
    if (jt>=Maxy-1) jt -= 1e-10;
    if (it<=0) it += 1e-10;
    if (it>=Maxx-1) it -= 1e-10;
    // get integer parts:
    long int jti = (long int)floor(jt);
    long int iti = (long int)floor(it);

// Cubic interpolation
    if (jti<0) jti=0;
    if (jti>=Maxy-4) 
      {
        jti=Maxy-4; // need 4 points
        // cout<<"out of y boundary"<<endl;
      }

    if (iti<0) iti=0;       
    if (iti>=Maxx-4) 
      { 
        iti=Maxx-4; // need 4 points
       // cout<<"out of x boundary"<<endl;
      }
//interpolation is done on the x-y grid, 

    double xfraction = it-iti;
    double yfraction = jt-jti;
     // row interpolation + extrapolation first
    double A0 = interpCubic4Points(unshiftedTable[iRap][iti][jti], unshiftedTable[iRap][iti][jti+1], 
      unshiftedTable[iRap][iti][jti+2], unshiftedTable[iRap][iti][jti+3], 1, yfraction);

    double A1 = interpCubic4Points(unshiftedTable[iRap][iti+1][jti], unshiftedTable[iRap][iti+1][jti+1], 
      unshiftedTable[iRap][iti+1][jti+2], unshiftedTable[iRap][iti+1][jti+3], 1, yfraction);

    double A2 = interpCubic4Points(unshiftedTable[iRap][iti+2][jti], unshiftedTable[iRap][iti+2][jti+1], 
      unshiftedTable[iRap][iti+2][jti+2], unshiftedTable[iRap][iti+2][jti+3], 1, yfraction);

    double A3 = interpCubic4Points(unshiftedTable[iRap][iti+3][jti], unshiftedTable[iRap][iti+3][jti+1], 
      unshiftedTable[iRap][iti+3][jti+2], unshiftedTable[iRap][iti+3][jti+3], 1, yfraction);
    
    return interpCubic4Points(A0,A1,A2,A3,1, xfraction,true);  //set this to true to avoid 
                                                                //interpolate to negative value

    //revised in Mar.18, 2013
    //switch to linear interpolation since:
    //1. cubic interpolation assumes smooth function, which is unknown for fKLN output
    //2. negative result is given by cubic method
    //interp on boundary
  //       if (jti<0) jti=0;
  //       if (jti>=Maxy-2) 
  //         {
  //           jti=Maxy-2; // need 2 points
  //           // cout<<"out of y boundary"<<endl;
  //         }

  //       if (iti<0) iti=0;       
  //       if (iti>=Maxx-2) 
  //         { 
  //           iti=Maxx-2; // need 2 points
  //          // cout<<"out of x boundary"<<endl;
  //         }
 	// return Bilinear2dInterp(iti, jti, 1, 1, unshiftedTable[iRap][iti][jti], 
  //                  unshiftedTable[iRap][iti][jti+1], unshiftedTable[iRap][iti+1][jti+1], 
  //                  unshiftedTable[iRap][iti+1][jti]);
    }
}




void FreeStrm::ShiftDensity(const int iRap, double phip)
{
/* get the new profile after it shifts to one angle
*/
  double Phip=phip;

  //store intermediate table
  shiftedTable  = new double** [nRap];
  for(int iy=0;iy<nRap;iy++) {
  shiftedTable[iy] =  new double* [Maxx];
  for(int i=0;i<Maxx;i++) {
      shiftedTable[iy][i] = new double [Maxy];
      for(int j=0;j<Maxy;j++) {
          shiftedTable[iy][i][j]=0.0;
      }
    }
  }

  for(int i=0;i<Maxx;i++)
    for(int j=0;j<Maxy;j++)
    {
      shiftedTable[iRap][i][j]=GetDensity(iRap,i,j,Phip);
    }              
//  cout<<"Profile has been shifted to: "<<Phip<<endl;
}



void FreeStrm::OutputTable(const char *filename, const int iRap)
{
  ofstream of;
  of.open(filename, std::ios_base::out);

  double rap = rapMin+(rapMax-rapMin)/nRap*iRap;


//  Output dNd2rdy Table
    for(int i=0;i<Maxx;i++)
      for(int j=0;j<Maxy;j++)
    {
      of <<  setprecision(3) << setw(10) <<  rap 
          << setprecision(3) << setw(10) <<  Xmin+i*dx 
          << setprecision(3) << setw(10) <<  Ymin+j*dy
          << setprecision(12) << setw(22) << shiftedTable[iRap][i][j]
          << endl;
    }
  cout<<"Free Streamed Profile complete!"<<endl;
  of.close();
}




void FreeStrm::dumpBlockTable(const char *filename, double ***data, const int iRap)
{
  
  ofstream of;
  of.open(filename, std::ios_base::out);

    for(int i=0;i<Maxx;i++)
    {
       for(int j=0;j<Maxy;j++)
           {
               of <<scientific << setw(22) <<data[iRap][i][j];              
           }
    of << endl;
    }

  cout<<"Block Table complete!"<<endl;
  of.close();

}

//Following functions are used for debugging and testing if free streaming
//generate the desirable profile.
void FreeStrm::CreateDataTable(const char *filename, const int iRap)
{
  ofstream of;
  of.open(filename, std::ios_base::out);
  cout << "Work in test mode: using Gaussian Profile!" << endl;

  double rap = rapMin+(rapMax-rapMin)/nRap*iRap;

  for(int i=0;i<Maxx;i++)
  for(int j=0;j<Maxy;j++)
  for(int ipt=0;ipt<MaxPT;ipt++)
  {
    double x=Xmin+i*dx;
    double y=Ymin+j*dy;
    double ptstep = PTmin + ipt*dpt;
    of <<  setprecision(3) << setw(10) <<  rap 
        << setprecision(3) << setw(10) <<  x 
        << setprecision(3) << setw(10) <<  y 
        << setprecision(3) << setw(10) <<  ptstep 
        << setprecision(12) << setw(22) << GaussProfile(iRap,i,j,ipt)   //debugging
        << endl;
  }
  cout<<"Gaussian initial condition is created!"<<endl;
  of.close();
}


double FreeStrm::GaussProfile(int iRap, int i, int j, int ipt)
{
    double x=Xmin+i*dx;
    double y=Ymin+j*dy;
    double ptstep = PTmin + ipt*dpt;
    double rap_factor = iRap;

    rap_factor = 1.; //rapidity =0 now
    double fxy;
    // prefactor=1/(12.087747372898045e0), Gaussian normalization factor

    fxy= exp(-(x - 0.1)*(x - 0.1)/3.0 - (y - 0.1)*(y - 0.1)/3.0 - ptstep*ptstep)/12.087747372898045;


    return fxy;
}

double FreeStrm::BoxProfile(const int iRap, int i, int j, int ipt)
{   
    double x=Xmin+i*dx;
    double y=Ymin+j*dy;
    double ptstep = PTmin + ipt*dpt;
    double fxy;
    double prefactor=1/18.0;

    double rap_factor = iRap;
    rap_factor = 1.; //rapidity =0 now
    
    fxy= prefactor*(stepfunc(x+3)-stepfunc(x-3))*(stepfunc(y+3)-stepfunc(y-3))*exp(-ptstep*ptstep);

    return fxy;
}

  
double FreeStrm::getEpx(int nth_order, const int iRap)
// Epx = \int{dx*dy*e(x,y)*(y^2-x^2)*gamma(ux,uy)}/\int{dx*dy*e(x,y)*(y^2+x^2)*gamma(ux,uy)}
// gamma factor is inlucded, since Landau matching generates flow, in order to transform energy
// density to the lab frame, gamma factor should be included.
{
  cout<<"Calculating spatial eccentricity Epx----------------"<<endl;
  double Epx=0.0;

  double epx_nu_real = 0.0;  //numerator of epx
  double epx_nu_img = 0.0;
  double Epx_angle = 0.0;
  double epx_dn = 0.0;  //denominator of epx
  
  for(int i=0; i<Maxx; i++)
  {
      for(int j=0; j<Maxy; j++)
     {
       double x = Xmin + dx*i - Xcm;  //center the profile
       double y = Ymin + dy*j - Ycm;
       double r = sqrt(x*x + y*y);
       double phi = atan2(y,x);

       double ed_temp = shiftedTable[iRap][i][j];

       epx_nu_real += pow(r, nth_order) 
                  * cos(nth_order * phi) * ed_temp * dx*dy;
       epx_nu_img += pow(r, nth_order) 
                  * sin(nth_order * phi) * ed_temp * dx*dy;
       epx_dn += ed_temp * pow(r, nth_order)* dx*dy;
     }       
  }//<-> for i=0:Maxx
  
  Epx = sqrt( epx_nu_real * epx_nu_real + epx_nu_img * epx_nu_img )
       /(epx_dn + 1e-18);
  Epx_angle = -atan2( epx_nu_img, (epx_nu_real + 1e-18))/nth_order;
/*
  cout<<"Spatial Eccentricity complete!"<<endl;
  cout << "Epx_angle =" <<Epx_angle << endl;*/

  return Epx;
}
