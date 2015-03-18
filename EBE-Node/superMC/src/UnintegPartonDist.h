#ifndef UnintegPartonDist_h
#define UnintegPartonDist_h

class UnintegPartonDist
{
 private:
  int largeXflag;  // allow or not extrapolation of uGD to x>x0
 public:
  UnintegPartonDist() {largeXflag=1;}
  void setlgXflag(int flag) {largeXflag=flag;}
  int probelgXflag() {return largeXflag;}

  virtual ~UnintegPartonDist() {}
  virtual double getFunc(double qs2, double x, double pt2,double alp)=0;
  virtual double getFuncNF(double qs2, double x, double pt2)=0;
  virtual double getQs(double, double)=0;
  virtual double getMax_kt() =0;  // max transv. momentum for given UGD set

  static const double x0 = 0.01;  // assumed starting point of small-x evolution
  static const double xCut = 0.1;  // minimal x_proj for Large_x()
  static const double lgXlambda = 0.3; /* evolution speed of Qs(x), for large-x
				  extrapolation according to ~exp(lambda*Y) */
};
#endif
