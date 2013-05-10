#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<stdlib.h>
#include "Arsenal.h"

double Simpson_sum(double* array, int num, double h_step)
{
      double sum=0;
      for(int j=0;j<num-2;j+=2)
         sum += array[j]*1 + array[j+1]*4 + array[j+2]*1; 
      sum = sum*h_step/3.0;
      return(sum);
}

//**********************************************************************
vector< vector<double>* >* readBlockData(istream &stream_in)
// Return a nested vector of vector<double>* object. Each column of data
// is stored in a vector<double> array and the collection is the returned
// object. Data are read from the input stream "stream_in". Each line
// of data is processed by the stringToDoubles function. Note that the
// data block is dynamically allocated and is not release within the
// function.
// Note that all "vectors" are "new" so don't forget to delete them.
// Note also that the last line of data needs a return in order to be read.
{
  vector< vector<double>* >* data;
  vector<double> valuesInEachLine;
  long lineSize;
  long i; // temp variable
  char buffer[99999]; // each line should be shorter than this

  // first line:
  stream_in.getline(buffer,99999);
  valuesInEachLine = stringToDoubles(buffer);
  // see if it is empty:
  lineSize = valuesInEachLine.size();
  if (lineSize==0)
  {
    // empty:
    cout << "readBlockData warning: input stream has empty first row; no data read" << endl;
    return NULL;
  }
  else
  {
    // not empty; allocate memory:
    data = new vector< vector<double>* >(lineSize);
    for (i=0; i<lineSize; i++) (*data)[i] = new vector<double>;
  }

  // rest of the lines:
  while (stream_in.eof()==false)
  {
    // set values:
    for (i=0; i<lineSize; i++) (*(*data)[i]).push_back(valuesInEachLine[i]);
    // next line:
    stream_in.getline(buffer,99999);
    valuesInEachLine = stringToDoubles(buffer);
  }

  return data;
}

//**********************************************************************
void releaseBlockData(vector< vector<double>* >* data)
// Use to delete the data block allocated by readBlockData function.
{
  if (data)
  {
    for (unsigned long i=0; i<data->size(); i++) delete (*data)[i];
    delete data;
  }
}


//**********************************************************************
vector<double> stringToDoubles(string str)
// Return a vector of doubles from the string "str". "str" should
// be a string containing a line of data.
{
  stringstream sst(str+" "); // add a blank at the end so the last data will be read
  vector<double> valueList;
  double val;
  sst >> val;
  while (sst.eof()==false)
  {
    valueList.push_back(val);
    sst >> val;
  }
  return valueList;
}


//**********************************************************************
double interpCubicDirect(vector<double>* x, vector<double>* y, double x0)
// Returns the interpreted value of y=y(x) at x=x0 using cubic polynomial interpolation method.
// -- x,y: the independent and dependent double x0ables; x is assumed to be equal spaced and increasing
// -- x0: where the interpolation should be performed
{
  long size = x->size();
  if (size==1) {cout<<"interpCubicDirect warning: table size = 1"; return (*y)[0];}
  double dx = (*x)[1]-(*x)[0]; // increment in x

  // if close to left end:
  if (abs(x0-(*x)[0])<dx*1e-30) return (*y)[0];

  // find x's integer index
  long idx = floor((x0-(*x)[0])/dx);

  if (idx<0 || idx>=size-1)
  {
    cout    << "interpCubicDirect: x0 out of bounds." << endl
            << "x ranges from " << (*x)[0] << " to " << (*x)[size-1] << ", "
            << "x0=" << x0 << ", " << "dx=" << dx << ", " << "idx=" << idx << endl;
    exit(1);
  }

  if (idx==0)
  {
    // use quadratic interpolation at left end
    double A0 = (*y)[0], A1 = (*y)[1], A2 = (*y)[2], deltaX = x0 - (*x)[0]; // deltaX is the increment of x0 compared to the closest lattice point
    return (A0-2.0*A1+A2)/(2.0*dx*dx)*deltaX*deltaX - (3.0*A0-4.0*A1+A2)/(2.0*dx)*deltaX + A0;
  }
  else if (idx==size-2)
  {
    // use quadratic interpolation at right end
    double A0 = (*y)[size-3], A1 = (*y)[size-2], A2 = (*y)[size-1], deltaX = x0 - ((*x)[0] + (idx-1)*dx);
    return (A0-2.0*A1+A2)/(2.0*dx*dx)*deltaX*deltaX - (3.0*A0-4.0*A1+A2)/(2.0*dx)*deltaX + A0;
  }
  else
  {
    // use cubic interpolation
    double A0 = (*y)[idx-1], A1 = (*y)[idx], A2 = (*y)[idx+1], A3 = (*y)[idx+2], deltaX = x0 - ((*x)[0] + idx*dx);
    //cout << A0 << "  " << A1 << "  " << A2 << "  " << A3 << endl;
    return (-A0+3.0*A1-3.0*A2+A3)/(6.0*dx*dx*dx)*deltaX*deltaX*deltaX
            + (A0-2.0*A1+A2)/(2.0*dx*dx)*deltaX*deltaX
            - (2.0*A0+3.0*A1-6.0*A2+A3)/(6.0*dx)*deltaX
            + A1;
  }

}




//**********************************************************************
double interpLinearDirect(vector<double>* x, vector<double>* y, double x0)
// Returns the interpreted value of y=y(x) at x=x0 using linear interpolation method.
// -- x,y: the independent and dependent double x0ables; x is assumed to be equal spaced and increasing
// -- x0: where the interpolation should be performed
{
  long size = x->size();
  if (size==1) {cout<<"interpLinearDirect warning: table size = 1"<<endl; return (*y)[0];}
  double dx = (*x)[1]-(*x)[0]; // increment in x

  // if close to left end:
  if (abs(x0-(*x)[0])<dx*1e-30) return (*y)[0];

  // find x's integer index
  long idx = floor((x0-(*x)[0])/dx);

  if (idx<0 || idx>=size-1)
  {
    cout    << "interpLinearDirect: x0 out of bounds." << endl
            << "x ranges from " << (*x)[0] << " to " << (*x)[size-1] << ", "
            << "x0=" << x0 << ", " << "dx=" << dx << ", " << "idx=" << idx << endl;
    exit(1);
  }

  return (*y)[idx] + ((*y)[idx+1]-(*y)[idx])/dx*(x0-(*x)[idx]);

}




//**********************************************************************
double interpNearestDirect(vector<double>* x, vector<double>* y, double x0)
// Returns the interpreted value of y=y(x) at x=x0 using nearest interpolation method.
// -- x,y: the independent and dependent double x0ables; x is assumed to be equal spaced and increasing
// -- x0: where the interpolation should be performed
{
  long size = x->size();
  if (size==1) {cout<<"interpNearestDirect warning: table size = 1"<<endl; return (*y)[0];}
  double dx = (*x)[1]-(*x)[0]; // increment in x

  // if close to left end:
  if (abs(x0-(*x)[0])<dx*1e-30) return (*y)[0];

  // find x's integer index
  long idx = floor((x0-(*x)[0])/dx);

  if (idx<0 || idx>=size-1)
  {
    cout    << "interpNearestDirect: x0 out of bounds." << endl
            << "x ranges from " << (*x)[0] << " to " << (*x)[size-1] << ", "
            << "x0=" << x0 << ", " << "dx=" << dx << ", " << "idx=" << idx << endl;
    exit(1);
  }

  return x0-(*x)[idx]>dx/2 ? (*y)[idx+1] : (*y)[idx];

}




//**********************************************************************
double interpCubicMono(vector<double>* x, vector<double>* y, double xx)
// Returns the interpreted value of y=y(x) at x=x0 using cubic polynomial interpolation method.
// -- x,y: the independent and dependent double x0ables; x is *NOT* assumed to be equal spaced but it has to be increasing
// -- xx: where the interpolation should be performed
{
  long size = x->size();
  if (size==1) {cout<<"interpCubicMono warning: table size = 1"<<endl; return (*y)[0];}

  // if close to left end:
  if (abs(xx-(*x)[0])<((*x)[1]-(*x)[0])*1e-30) return (*y)[0];

  // find x's integer index
  long idx = binarySearch(x, xx);

  if (idx<0 || idx>=size-1)
  {
    cout    << "interpCubicMono: x0 out of bounds." << endl
            << "x ranges from " << (*x)[0] << " to " << (*x)[size-1] << ", "
            << "xx=" << xx << ", " << "idx=" << idx << endl;
    exit(1);
  }

  if (idx==0)
  {
    // use linear interpolation at the left end
    return (*y)[0] + ( (*y)[1]-(*y)[0] )/( (*x)[1]-(*x)[0] )*( xx-(*x)[0] );
  }
  else if (idx==size-2)
  {
    // use linear interpolation at the right end
    return (*y)[size-2] + ( (*y)[size-1]-(*y)[size-2] )/( (*x)[size-1]-(*x)[size-2] )*( xx-(*x)[size-2] );
  }
  else
  {
    // use cubic interpolation
    long double y0 = (*y)[idx-1], y1 = (*y)[idx], y2 = (*y)[idx+1], y3 = (*y)[idx+2];
    long double y01=y0-y1, y02=y0-y2, y03=y0-y3, y12=y1-y2, y13=y1-y3, y23=y2-y3;
    long double x0 = (*x)[idx-1], x1 = (*x)[idx], x2 = (*x)[idx+1], x3 = (*x)[idx+2];
    long double x01=x0-x1, x02=x0-x2, x03=x0-x3, x12=x1-x2, x13=x1-x3, x23=x2-x3;
    long double x0s=x0*x0, x1s=x1*x1, x2s=x2*x2, x3s=x3*x3;
    long double denominator = x01*x02*x12*x03*x13*x23;
    long double C0, C1, C2, C3;
    C0 = (x0*x02*x2*x03*x23*x3*y1
          + x1*x1s*(x0*x03*x3*y2+x2s*(-x3*y0+x0*y3)+x2*(x3s*y0-x0s*y3))
          + x1*(x0s*x03*x3s*y2+x2*x2s*(-x3s*y0+x0s*y3)+x2s*(x3*x3s*y0-x0*x0s*y3))
          + x1s*(x0*x3*(-x0s+x3s)*y2+x2*x2s*(x3*y0-x0*y3)+x2*(-x3*x3s*y0+x0*x0s*y3))
          )/denominator;
    C1 = (x0s*x03*x3s*y12
          + x2*x2s*(x3s*y01+x0s*y13)
          + x1s*(x3*x3s*y02+x0*x0s*y23-x2*x2s*y03)
          + x2s*(-x3*x3s*y01-x0*x0s*y13)
          + x1*x1s*(-x3s*y02+x2s*y03-x0s*y23)
          )/denominator;
    C2 = (-x0*x3*(x0s-x3s)*y12
          + x2*(x3*x3s*y01+x0*x0s*y13)
          + x1*x1s*(x3*y02+x0*y23-x2*y03)
          + x2*x2s*(-x3*y01-x0*y13)
          + x1*(-x3*x3s*y02+x2*x2s*y03-x0*x0s*y23)
          )/denominator;
    C3 = (x0*x03*x3*y12
          + x2s*(x3*y01+x0*y13)
          + x1*(x3s*y02+x0s*y23-x2s*y03)
          + x2*(-x3s*y01-x0s*y13)
          + x1s*(-x3*y02+x2*y03-x0*y23)
          )/denominator;
/*    cout  << x0s*x03*x3s*y12 << "  "
          <<  x2*x2s*(x3s*y01+x0s*y13) << "   "
          <<  x1s*(x3*x3s*y02+x0*x0s*y23-x2*x2s*y03) << "  "
          <<  x2s*(-x3*x3s*y01-x0*x0s*y13) << "  "
          <<  x1*x1s*(-x3s*y02+x2s*y03-x0s*y23) << endl;
    cout << denominator << endl;

    cout << x0 << " " << x1 << "  " << x2 << "  " << x3 << endl;
    cout << y0 << " " << y1 << "  " << y2 << "  " << y3 << endl;
    cout << C0 << "  " << C1 << "  " << C2 << "  " << C3 << endl;*/
    return C0 + C1*xx + C2*xx*xx + C3*xx*xx*xx;
  }

}




//**********************************************************************
double interpLinearMono(vector<double>* x, vector<double>* y, double xx)
// Returns the interpreted value of y=y(x) at x=x0 using linear interpolation method.
// -- x,y: the independent and dependent double x0ables; x is *NOT* assumed to be equal spaced but it has to be increasing
// -- xx: where the interpolation should be performed
{
  long size = x->size();
  if (size==1) {cout<<"interpLinearMono warning: table size = 1"<<endl; return (*y)[0];}

  // if close to left end:
  if (abs(xx-(*x)[0])<((*x)[1]-(*x)[0])*1e-30) return (*y)[0];

  // find x's integer index
  long idx = binarySearch(x, xx);

  if (idx<0 || idx>=size-1)
  {
    cout    << "interpLinearMono: x0 out of bounds." << endl
            << "x ranges from " << (*x)[0] << " to " << (*x)[size-1] << ", "
            << "xx=" << xx << ", " << "idx=" << idx << endl;
    exit(1);
  }

  return (*y)[idx] + ( (*y)[idx+1]-(*y)[idx] )/( (*x)[idx+1]-(*x)[idx] )*( xx-(*x)[idx] );

}




//**********************************************************************
double interpNearestMono(vector<double>* x, vector<double>* y, double xx)
// Returns the interpreted value of y=y(x) at x=x0 using nearest interpolation method.
// -- x,y: the independent and dependent double x0ables; x is *NOT* assumed to be equal spaced but it has to be increasing
// -- xx: where the interpolation should be performed
{
  long size = x->size();
  if (size==1) {cout<<"interpNearestMono warning: table size = 1"<<endl; return (*y)[0];}

  // if close to left end:
  if (abs(xx-(*x)[0])<((*x)[1]-(*x)[0])*1e-30) return (*y)[0];

  // find x's integer index
  long idx = binarySearch(x, xx);

  if (idx<0 || idx>=size-1)
  {
    cout    << "interpNearestMono: x0 out of bounds." << endl
            << "x ranges from " << (*x)[0] << " to " << (*x)[size-1] << ", "
            << "xx=" << xx << ", " << "idx=" << idx << endl;
    exit(1);
  }

  return xx-(*x)[idx] > (*x)[idx+1]-xx ? (*y)[idx+1] : (*y)[idx];

}




//**********************************************************************
double invertFunc(double (*func)(double), double y, double xL, double xR, double dx, double x0, double relative_accuracy)
//Purpose:
//  Return x=func^(-1)(y) using Newton method.
//  -- func: double 1-argument function to be inverted
//  -- xL: left boundary (for numeric derivative)
//  -- xR: right boundary (for numeric derivative)
//  -- dx: step (for numeric derivative)
//  -- x0: initial value
//  -- y: the value to be inverted
//  -- Returns inverted value
//Solve: f(x)=0 with f(x)=table(x)-y => f'(x)=table'(x)
{
  double accuracy;
  int tolerance;

  double XX1, XX2; // used in iterations
  double F0, F1, F2, F3, X1, X2; // intermedia variables
  int impatience; // number of iterations


  // initialize parameters
  accuracy = dx*relative_accuracy;

  tolerance = 60;
  impatience = 0;

  // initial value, left point and midxle point
  XX2 = x0;
  XX1 = XX2-10*accuracy; // this value 10*accuracy is meanless, just to make sure the check in the while statement goes through

  while (abs(XX2-XX1)>accuracy)
  {
    XX1 = XX2; // copy values

    // value of function at XX
    F0 = (*func)(XX1) - y; // the value of the function at this point

    // decide X1 and X2 for differentiation
    if (XX1>xL+dx)
      X1 = XX1 - dx;
    else
      X1 = xL;

    if (XX1<xR-dx)
      X2 = XX1 + dx;
    else
      X2 = xR;

    // get values at X1 and X2
    F1 = (*func)(X1);
    F2 = (*func)(X2);
    F3 = (F1-F2)/(X1-X2); // derivative at XX1

    XX2 = XX1 - F0/F3; // Newton's mysterious method

    impatience = impatience + 1;
    //cout << "impatience=" << impatience << endl;
    if (impatience>tolerance)
    {
      cout << "invertFunc: " << "max number of iterations reached." << endl;
      exit(-1);
    }

  } // <=> abs(XX2-XX1)>accuracy

  return XX2;
}

//**********************************************************************
vector<double> *zq_x_global, *zq_y_global;
double invertTableDirect_hook(double xx) {return interpCubicDirect(zq_x_global,zq_y_global,xx);}
double invertTableDirect(vector<double>* x, vector<double>* y, double y0, double x0, double relative_accuracy)
// Return x0=y^(-1)(y0) for y=y(x); use interpCubic and invertFunc.
//  -- x,y: the independent and dependent variables. x is assumed to be equal-spaced.
//  -- y0: where the inversion should be performed.
//  -- x0: initial guess
{
  long size = x->size();
  if (size==1) return (*y)[0];
  zq_x_global = x; zq_y_global = y;
  return invertFunc(&invertTableDirect_hook, y0, (*x)[0], (*x)[size-1], (*x)[1]-(*x)[0], x0, relative_accuracy);
}

//**********************************************************************
long binarySearch(vector<double>* A, double value)
// Return the index of the largest number less than value in the list A
// using binary search.
{
   int length = A->size();
   int idx_i, idx_f, idx;
   idx_i = 0;
   idx_f = length-1;
   if(value > (*A)[idx_f])
   {
      cout << "binarySearch: desired value is too large,  exceed the end of the table." << endl;
      exit(1);
   }
   if(value < (*A)[idx_i])
   {
      cout << "binarySearch: desired value is too small, exceed the begin of table." << endl;
      exit(1);
   }
   idx = (int) floor((idx_f+idx_i)/2.);
   while((idx_f-idx_i) > 1)
   {
     if((*A)[idx] < value)
        idx_i = idx;
     else
        idx_f = idx;
     idx = (int) floor((idx_f+idx_i)/2.);
   }
   return(idx_i);
}

void outputFunctionerror(string function_name, string massage, double value, int level)
{
   switch(level)
   {
      case 0:
         break;
      case 1:
         cout << ">>> Warning: " << function_name << endl; 
         cout << ">>> Message: " << massage << endl;
         cout << ">>> Related value: " << value << endl;
         break;
      case 2:
         cout << ">>> Error: " << function_name << endl; 
         cout << ">>> Message: " << massage << endl;
         cout << ">>> Related value: " << value << endl;
         break;
      case 3:
         cout << ">>> Fatal Error: " << function_name << endl; 
         cout << ">>> Message: " << massage << endl;
         cout << ">>> Related value: " << value << endl;
         exit(1);
         break;
      default:
         cout << ">>> Error: outputFunctionerror:" << endl;
         cout << ">>> Message: wrong errorLevel" << endl;
         cout << ">>> Related value: " << level << endl;
         exit(1);
         break;
   }
}

//**********************************************************************
string toLower(string str)
// Convert all character in string to lower case
{
  string tmp = str;
  for (string::iterator it=tmp.begin(); it<=tmp.end(); it++) *it = tolower(*it);
  return tmp;
}
//**********************************************************************
string trim(string str)
// Convert all character in string to lower case
{
  string tmp = str;
  long number_of_char = 0;
  for (size_t ii=0; ii<str.size(); ii++)
    if (str[ii]!=' ' && str[ii]!='\t')
    {
      tmp[number_of_char]=str[ii];
      number_of_char++;
    }
  tmp.resize(number_of_char);
  return tmp;
}

//**********************************************************************
double stringToDouble(string str)
// Return the 1st doubles number read from the string "str". "str" should be a string containing a line of data.
{
  stringstream sst(str+" "); // add a blank at the end so the last data will be read
  double val;
  sst >> val;
  return val;
}

