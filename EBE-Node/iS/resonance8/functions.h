/* functions.h        header file for function.c
*  Sollfrank Josef           Nov .98                      */

#ifndef function_h
#define function_h

// SB code
#define	NUMDECAY	1700 /* size of array for storage of the decays */
#define	NUMPARTICLE	400 /*  size of array for storage of the particles */
#define MAXINTV         20000000 /* size of arry for Montecarlo numbers */
#define MHALF           (MAXINTV/2)
#define NPT		50 /* size of arry for storage of the pt-spectrum */
#define NPHI		120 /* size of arry for storage of the y-spectrum */
#define NPHI1		NPHI + 1
#define	PI		3.14159265358979323 /* any question? */
#define ALPHA           0.00729735308 /* fine structure constant 1/137.035...*/
#define HBARC		0.197327054 /* = plank constant times speed of light */
#define	HBARC3		(HBARC*HBARC*HBARC)
#define FILEDIM         140

void   readin(char filename[FILEDIM], int* particlemax, int* decaymax);
void   readspec(char  specfile[FILEDIM], int *particlemax, int *decaymax);
void   writespec(int particlemax, char outdir[FILEDIM]);
double Edndp3(double yr, double ptr, double phirin, int res_num);

#endif
