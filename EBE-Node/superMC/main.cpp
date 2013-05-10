#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <cstring>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <ctime>
#include <unistd.h>
#include <sys/time.h>
#include "MakeDensity.h"
#include "ParamDefs.h"
#include "ParameterReader.h"
using namespace std;

//###########################################################################

int main(int argc, char *argv[]) {

  // Read-in parameters
  ParameterReader paraRdr;
  paraRdr.readFromFile("parameters.dat");
  paraRdr.readFromArguments(argc, argv);
  paraRdr.echo();

  // init random seed from system time
  timeval a;
  gettimeofday(&a, 0);
  int randomSeed=paraRdr.getVal("randomSeed");
  if (randomSeed<0) randomSeed=a.tv_usec; // randomSeed<0 means to use CPU clock
  srand48(randomSeed);

  MakeDensity *dens = new MakeDensity(&paraRdr);

  int nevent = paraRdr.getVal("nev");       // # of events
  int operation = paraRdr.getVal("operation");

  time_t start, end;
  double cpu_time_used;
  start = clock();

  switch(operation)
  {
    case 1:
      // generate a set of profiles, used for e-by-e calculations
      dens->generate_profile_ebe(nevent);
      break;
    case 2:
      // generate a set of profiles, used for e-by-e calculations with Jet quenching
      dens->generate_profile_ebe_Jet(nevent);
      break;
    case 3:
      // generate a set of averaged profiles, for entropy and/or energy densities (controallable using parameters.dat)
      dens->generate_profile_average(nevent);
      break;
    case 9:
      // generate a table of eccentricities
      dens->generateEccTable(nevent);
      break;
    default:
      cout << "Error: operation choice " << operation << " not recognized." << endl;
  }
  end = clock();
  cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;

  cout << "Time elapsed (in seconds): " << cpu_time_used << endl;

  if (dens) delete dens;

  return 0;
}
