#include <cmath>
#include <iostream>
using namespace std;

namespace hadronxsec {

//**********************************************************************
//...fit for high energy cross sections by regge theory
//...(particle data group  1996 Phys. Rev. D54)
//...sig_tot=X*s^eta + Y*s^eps
//   0: pp    total    (50GeV/c)
//   1: p~ p  total    (50GeV/c)
//   2: np    total    (50GeV/c)
//   3: p~ n  total    (50GeV/c)
//   4: pd    total    (50GeV/c)
//   5: p~ d total     (50GeV/c)
//   6: pi+ p total    (10GeV/c)
//   7: pi- p total    (10GeV/c)
//   8: pi+- d total   (10GeV/c)
//   9: k+ p  total    (10GeV/c)
//  10: k+ n  total    (10GeV/c)
//  11: k- p  total    (10GeV/c)
//  12: k- n  total    (10GeV/c)
//  13: k+ d  totol    (10GeV/c)
//  14: k- d  total    (10GeV/c)
//  15: gamm p total   (12GeV/c)
double totalXsection(const double srt,const int i)
{
    static const double eta[16]={0.46,0.46,0.46,0.46,
      0.45,0.45, 0.45,0.45,  0.43, 0.5,0.5,0.5,0.5,    0.47,0.47, 0.46};
    static const double eps[16]={0.079,0.079,0.079,0.079,
	0.09,0.09,0.079,0.079,0.088,
        0.079,0.079,0.079,0.079, 0.082,0.082,0.075};
    static const double x[16]={22.0,22.0, 22.3,22.3, 35.7,35.7, 13.7,13.7, 23.2,
	12.2,12.2,12.2,12.2,  21.7,21.7, 0.071};
    static const double y[16]={56.1, 98.2, 55.0, 92.7,  179.0, 270.6,
       27.8,35.9, 85.5, 8.3,8.3,  26.4,26.4, 26.2,64.8, 0.12};

      return x[i]*pow(srt*srt,eps[i]) + y[i]*pow(srt*srt,-eta[i]);
}

double elasticXsection(double sig, double srt, int i1, int i2)
{
    const double bhad[3] = {2.3,1.4,0.23}; //1:baryon,2:meson,3:J/psi
    const double eps=0.0808;
    const double facel=0.0511;

    double bel=2.0*bhad[i1]+2.0*bhad[i2] + 4.0*pow(srt*srt,eps) - 4.2;
    return facel*sig*sig/bel;
}

#ifdef MAIN
int main()
{
    double srt = 14000;
    double sig = totalHadronicXsection(srt,0);
    double sigel = elasticHadronicXsection(sig,srt,0,0);

    cout << sig << " sigel= " << sigel << endl;
    cout << sig-sigel << endl;

    return 0;
}
#endif

}  // namespae hadronxsec
