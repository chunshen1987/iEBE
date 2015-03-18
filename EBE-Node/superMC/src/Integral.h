#ifndef Integral_h
#define Integral_h

class Integral
{
private:
    virtual double evalPar(const double x)=0;
public:
    double Gauss(double a, double b, double epsilon);
};

#endif // Integral_h
