#ifndef ARSENAL_H
#define ARSENAL_H

#include<fstream>
#include<string>
#include<vector>

using namespace std;

double Simpson_sum(double* , int, double);

vector< vector<double>* >* readBlockData(istream &stream_in);
void releaseBlockData(vector< vector<double>* >* data);
vector<double> stringToDoubles(string str);

double interpCubicDirect(vector<double>* x, vector<double>* y, double xx);
double interpCubicMono(vector<double>* x, vector<double>* y, double xx);
double interpLinearDirect(vector<double>* x, vector<double>* y, double xx);
double interpLinearMono(vector<double>* x, vector<double>* y, double xx);
double interpNearestDirect(vector<double>* x, vector<double>* y, double xx);
double interpNearestMono(vector<double>* x, vector<double>* y, double xx);

double invertFunc(double (*func)(double), double y, double xL, double xR, double dx, double x0, double relative_accuracy=1e-10);

double invertTableDirect_hook(double xx);
double invertTableDirect(vector<double>* x, vector<double>* y, double y0, double x0, double relative_accuracy=1e-10);

long binarySearch(vector<double>* A, double value);

void outputFunctionerror(string function_name, string massage, double value, int level);

string toLower(string str);
string trim(string str);
double stringToDouble(string str);

#endif
