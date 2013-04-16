#ifndef int_h
#define int_h

double gauss(int n, double (*func)(double, void*), double xlo, double xhi, void *optvec);
/*
void gausspts( int n, double xlo, double xhi, double *x, double *w );
double	gaussn(int n, ndiv, double (*func)(), double xlo, double xhi, double *para );
double	gala(int n, double (*func)(), double xlo, double invslope, void *optvec);
double gahe(int n, double (*func)(), double center, double width, void *optvec);
double	gauche(int n,double (*func)(), double base, double pole, void *optvec );
double	gaussp();
double	galap();
double	gahep();

double	gaussnbyn();
*/

#endif