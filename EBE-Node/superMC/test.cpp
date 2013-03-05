#include <iostream>
#include <fstream>
#include <iomanip>
#include "Table.h"
#include "RandomVariable.h"
#include "NBD.h"
#include "arsenal.h"
#include "TableFunction.h"
#include "Stopwatch.h"

using namespace std;

int main()
{
  Stopwatch sw;
/*
  Table A("sampleTable.dat");
  cout << setprecision(12) << A.interpCubicFlat(1,2,4.6) << endl;
  cout << A.interpCubicMono(2,1,6.1) << endl;
  cout << A.invertMono(1,2,6.1) << endl;
  Table B("sampleTable1.dat");
  cout << B.invertMono(1,2,24.1) << endl;
  cout << A.invertMono(1,2,6.1) << endl;
  cout << A.getFirst(1) << "," << B.getLast(1) << endl;
*/
  Table A;
  //A.extendTable(3,4);
  A.set(3,4,0.3);
  formatedPrint(2, A.get(1,4), A.get(3,4));

  TableFunction f,g;
  f.setMappingTable(1,1,1);
  f.setMappingTable(2,2,4);
  f.setMappingTable(3,3,9);

  cout << f.map(2.5) << endl;


  RandomVariable randVar;
  randVar.pdfTab->loadMappingTableFromFile("gaussian.dat");
  ofstream of1("dat1.dat");
  randVar.constructEnvelopTab(0,0.2,9,9);
  for (double x=-1.5; x<1.5; x+=0.05)
    of1 << x << "  " << randVar.pdf(x) << "  " << randVar.envelopPdf(x) << endl;

  randVar.envelopPdfTab->printFunction();

  cout << randVar.calculateMoments(1,-1,1) << endl;
  cout << randVar.calculateMoments(2,-1,1) << endl;
  cout << randVar.calculateMoments(3,-1,1) << endl;



  NBD nbd(0.3,0.7);

  cout << nbd.envelopPdfTab->getXMin() << "," << nbd.envelopPdfTab->getXMax() << endl;
  cout << nbd.envelopInvCDFTab->getXMin() << "," << nbd.envelopInvCDFTab->getXMax() << endl;
  cout << nbd.envelopInvCDFTab->getYMin() << "," << nbd.envelopInvCDFTab->getYMax() << endl;

  cout << "yes" << endl;
  ofstream fs1("rand1.dat");
  sw.tic();
  for (double x=-0.1; x<=4; x+=0.02)
    fs1 << x << "  " << nbd.pdf(x) << "  "
        << nbd.envelopPdf(x) << endl;
  sw.toc();
  cout << "NBD sampling takes " << sw.takeTime() << "seconds." << endl;

  ofstream fs2("rand2.dat");
  sw.tic();
  for (int i=0; i<=100000; i++) fs2 << nbd.rand(0.3,0.7) << endl;
  sw.toc();
  cout << "NBD sampling takes " << sw.takeTime() << "seconds." << endl;

  /*
  sw.tic();
  ofstream fs1("rand1.dat");
  for (int i=0; i<=100000; i++) fs1 << randVar.sampleAccToPDFUsingInvCDF(20,-2,2) << endl;
  sw.toc();
  cout << "Method 1 takes " << sw.takeTime() << "seconds." << endl;

  sw.tic();
  ofstream fs2("rand2.dat");
  for (int i=0; i<=100000; i++) fs2 << randVar.sampleAccToPDFDirect(0.06,-2,2) << endl;
  sw.toc();
  cout << "Method 1 takes " << sw.takeTime() << "seconds." << endl;

  sw.tic();
  ofstream fs3("rand3.dat");
  for (int i=0; i<=100000; i++) fs1 << randVar.sampleAccToPDFUsingInvCDF(20,-2,2) << endl;
  sw.toc();
  cout << "Method 1 takes " << sw.takeTime() << "seconds." << endl;

  sw.tic();
  ofstream fs4("rand4.dat");
  for (int i=0; i<=100000; i++) fs2 << randVar.sampleAccToPDFDirect(0.06,-2,2) << endl;
  sw.toc();
  cout << "Method 1 takes " << sw.takeTime() << "seconds." << endl;
  */
}
