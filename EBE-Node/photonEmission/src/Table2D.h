#ifndef Table2D_h
#define Table2D_h

#include <string>
#include <vector>
using namespace std;

class Table2D
{
  private:
    std::vector< std::vector<double>* >* data;
    int tb_sizeX, tb_sizeY;
    double varX_min;
    double varY_min;
    double dvarX;
    double dvarY;

  public:
    Table2D();
    Table2D(string );
    void loadTableFromFile(string );
    void setVarXmin(double xmin) {varX_min = xmin;}
    void setVarYmin(double ymin) {varY_min = ymin;}
    void setdvarX(double dx) {dvarX = dx;}
    void setdvarY(double dy) {dvarY = dy;}
    double getVarXmin() {return(varX_min);}
    double getVarYmin() {return(varY_min);}
    double getdvarX() {return(dvarX);}
    double getdvarY() {return(dvarY);}
    int getTbsizeX() {return(tb_sizeX);}
    int getTbsizeY() {return(tb_sizeY);}
    double getTbdata(int idx, int idy) {return((*(*data)[idy])[idx]);}
    void setTbdata(int idx, int idy, double value) {(*(*data)[idy])[idx]=value;}
    void outputTabletoFile(string filename);
};

#endif

