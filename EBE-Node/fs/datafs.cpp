#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring>

using namespace std;
int main()
{	
//change the form of output file from 201*201 matrix to 40401 column;

	ifstream DataFile("matrixdata.dat");

    double result[40401];
        for(int i=0;i<40401;i++)
            	result[i]=0.0;

    cout<<"table created"<<endl;

	ofstream of;
    of.open("coldata.dat", std::ios_base::out);

    cout<<"table open"<<endl;


     while (!DataFile.eof())
     {
     	            for(int i=0;i<40401;i++)
                    	DataFile>>result[i];
                    	

    }
            for(int i=0;i<40401;i++)
    			of << scientific <<result[i]<<endl;
    cout<<"data read"<<endl;


	of.close();
	DataFile.close();


	cout<<"table converted!"<<endl;

	return 0;



}