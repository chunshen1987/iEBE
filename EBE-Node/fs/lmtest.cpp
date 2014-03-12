#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include <ctime>
#include <stdlib.h>
#include "LdMatching.h"
#include "Freestreaming.h"
#include "gauss_quadrature.h"
#include "ParameterReader.h"

using namespace std;

int main(int argc, char *argv[])
{
  ParameterReader Params;
  Params.readFromFile("LMParameters.dat");
  Params.readFromArguments(argc, argv);
  Params.echo();
  //switch between run mode: one event or multiple events
  int event_flag = Params.getVal("event_mode");
  int nevents = Params.getVal("nevents");
  int event_num;
  if(event_flag <= 0)
  {
    event_num = 1;  //run multiple eventss
  }
  else if(event_flag > 0)
  {
    event_num = event_flag;
    nevents = event_flag;   //only run one events
  }

  //Timing the current run
  time_t start, end;
  double cpu_time_used;
  start = clock();

  //processing events
  for( ;event_num<=nevents;event_num++)
  {
    //prepare readin filename for event-by-event eccentricity fluctuation
    ostringstream filename_stream;
    filename_stream.str("");
    filename_stream << "data/events/sd_event_"
                    << event_num  <<"_block.dat";
//    filename_stream << "data/sdAvg_order_2_block"
//                    << ".dat";

    //prepare data directory for final profiles of different events and 
    //different matching time
    ostringstream result_dir_stream;
    result_dir_stream.str("");
    result_dir_stream << "./data/result/event_" << event_num;
    string result_directory = result_dir_stream.str();
    system(("rm -rf " + result_directory).c_str());
    system(("mkdir " + result_directory).c_str());

    LdMatching *Matching;
    Matching = new LdMatching(&Params, result_directory);
    Matching->MultiMatching(filename_stream.str());  //do the matching
    delete Matching;
  }

  end = clock();
  cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;

  cout << "Time elapsed (in seconds): " << cpu_time_used << endl;
  return 0;
}
