/***********************************************************************
ArraySaver class is used to store or load an array onto a disk. The purpose is to use this class to backup DOUBLE type data so that if a program is terminated unexpectively, it can be restored to some extent.
An instance of the class should be associated to a file. If that file exists, the constructor will load data from that file; otherwise it will save to that file whenever the "snapshot" command is called. The file is deleted when the class is deleted (i.e. on normal exit). The snapshot can be manually loaded (not recommeded) by loadSnapshot method.
The grammar for initializing the class is bit complicated, see ArraySaver.cpp for details.
Version 1.1 (09-21-2011) Zhi Qiu
***********************************************************************/
#ifndef ArraySaver_header
#define ArraySaver_header

#include <fstream>
#include <string>
#include <cstdarg>
#include <iostream>

using namespace std;

template <class T>
class ArraySaver
{
  private:
    std::string filename;
    void* data;
    long depth; long* sizeList;
    void iterativeWrite(void* p, long level, std::ofstream &fp);
    void iterativeRead(void* p, long level, std::ifstream &fp);
  public:
    ArraySaver(std::string filename_in, void* data_in, int depth_in, ...);
    ~ArraySaver();
    void snapshot();
    void loadSnapshot(bool no_file_warning = true);
};



template <class type>
ArraySaver<type>::ArraySaver(string filename_in, void* data_in, int depth_in, ...)
/*
 -- filename_in: string of the associated filename.
 -- data_in: pointer to the data to be saved/loaded.
 -- depth_in: the "depth" of the pointer data_in. Note that the depth is not the dimension of the array, but the "reference level". For example double A[5][5] has depth 1 while double **A has depth 2.
 -- ...: list of size of each dimension of the data, the lower level goes first (to the left). For example, double **A = new double*[5], then 5 is the highest dimension, it should be put to the right.
*/
{
  filename = filename_in;
  data = data_in;
  depth = depth_in;
  sizeList = new long[depth];
  va_list ap;
  va_start(ap, depth_in); //Requires the last fixed parameter (to get the address)
  for(int j=0; j<depth; j++)
  {
    sizeList[j] = va_arg(ap, long); //Requires the type to cast to. Increments ap to the next argument.
  }
  va_end(ap);
  loadSnapshot(false);
}

template <class type>
ArraySaver<type>::~ArraySaver()
{
  delete[] sizeList;
  remove(filename.c_str());
};

template <class type>
void ArraySaver<type>::iterativeWrite(void* p, long level, ofstream &fp)
{
  if (level==1)
  {

    for (long i=0; i<sizeList[level-1]; i++)
    {
      fp << *((type*)p+i) << endl;
    }
  }
  else
  {
    for (long i=0; i<sizeList[level-1]; i++)
    {
      iterativeWrite(*((type**)p+i), level-1, fp);
    }
  }
}

template <class type>
void ArraySaver<type>::snapshot()
// Save all data to file with name "filename"
{
  ofstream fp(filename.c_str());
  iterativeWrite(data, depth, fp);
  fp.close();
}

template <class type>
void ArraySaver<type>::iterativeRead(void* p, long level, ifstream &fp)
{
  if (level==1)
  {

    for (long i=0; i<sizeList[level-1]; i++)
    {
      fp >> *((type*)p+i);
    }
  }
  else
  {
    for (long i=0; i<sizeList[level-1]; i++)
    {
      iterativeRead(*((type**)p+i), level-1, fp);
    }
  }
}

template <class type>
void ArraySaver<type>::loadSnapshot(bool no_file_warning)
// Load all data from file with name "filename"
// -- no_file_warning: if it is set to true (default), then an warning message is printed if the file does not exist; in both case the loading is skipped in such situations.
{
  ifstream fp(filename.c_str());
  if (fp) // load only if file open successfully
  {
    iterativeRead(data, depth, fp);
  }
  else if (no_file_warning)
  {
    cout << "ArraySaver::loadSnapshot warning: file " << filename << " does not exist." << endl;
  }
  fp.close();
}


#endif

/*----------------------------------------------------------------------
Change log:
09-21-2011: Ver 1.1:
  -- Due to the limitation of g++, the requirement of usage of templates forecs the combination of ArraySaver.cpp and ArraySaver.h file. From this version on, only the ArraySaver.h file needs to be included in the compiling process.

----------------------------------------------------------------------*/
