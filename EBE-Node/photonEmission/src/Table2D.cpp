#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include "stdlib.h"
#include <cmath>
#include "Arsenal.h"
#include "Table2D.h"


using namespace std;

Table2D::Table2D() {};

Table2D::Table2D(string filename) { loadTableFromFile(filename); };

//----------------------------------------------------------------------
void Table2D::loadTableFromFile(string data_filename)
// The Table data file (data_filename) is assumed to be a n-column file:
{
  ostringstream filename_stream;
  filename_stream << data_filename;
  fstream fs(filename_stream.str().c_str());
  if (fs.is_open()==false)
  {
    cout << "Table2D::loadTableFromFile error: the data file cannot be opened." << endl;
    cout << "Filename: " << data_filename << endl;
    exit(-1);
  }
  data = readBlockData(fs);
  tb_sizeY = (*data).size();
  tb_sizeX = (*(*data)[0]).size();
}

void Table2D::outputTabletoFile(string filename)
{
   ostringstream filename_stream;
   filename_stream << filename << ".dat";
   ofstream output(filename_stream.str().c_str());
   for(int i=0; i<tb_sizeX; i++)
   {
     for(int j=0; j<tb_sizeY; j++)
        output << scientific << setw(16) << setprecision(6)
               << (*(*data)[j])[i] << "   ";
     output << endl;
   }
}
