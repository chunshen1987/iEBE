/*	decay.c
	resonance decays for routines
	started 14 May 91, es and Jockel
*/
// Commented and adjusted by Evan Frodermann (Aug 2005)

#include	<string.h>
#include	<stdio.h>
#include	<math.h>
#include        <stdlib.h>
//#include	"reso.h"
#include	"decay.h"
#include	"int.h"


#define	PTS3	12		/* normalization of 3-body decay */
#define	PTS4	12		/* inv. mass integral 3-body    */

#define	PTN1	8
#define	PTN2	8		/* 2-body between the poles     */



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

extern struct par  {
	long int 	monval;     /* Montecarlo number according PDG */
	char	name[26];
	double	mass;
	double	width;
	int	gspin;      /* spin degeneracy */
	int	baryon;
	int	strange;
	int	charm;
	int	bottom;
	int	gisospin;  /* isospin degeneracy */
	int	charge;
	int	decays;    /* amount of decays listed for this resonance */
	int	stable;     /* defines whether this particle is considered
				as stable */
        int     npt;
        int     nphi;
	double  dNdptdphi[NPT][NPHI];
	double	mt[NPT];	/* mt values for spectrum */
	double	pt[NPT];	/* pt values for spectrum */
	double  slope[NPHI1];		/* assymtotic slope of mt-spectrum */
	} particle[NUMPARTICLE];

extern double  PHI[NPHI];    /* Array for phi-storage */

/* decay array for each decay listed */
extern struct de  {
	int	reso;       /* Montecarlo number of decaying resonance */
	int	numpart;    /* number of daughter particles after decay */
	double	branch;     /* branching ratio */
	int	part[5];    /* array of daughter particels Montecarlo number */
	}  particleDecay[NUMDECAY];

/* array for converting Montecarlo numbers in internal numbering of the
resonances */
extern int	partid[MAXINTV];

typedef struct pblockN
{
  double pt, mt, y, e, pl;	/* pt, mt, y of decay product 1 */
  double phi;
  double m1, m2, m3;		/* masses of decay products     */
  double mr;			/* mass of resonance            */
  double costh, sinth;
  double e0, p0;
  int res_num;			/* Montecarlo number of the Res. */
} pblockN;


typedef struct nblock
{
  double a, b, c, d;
} nblock;			/* for normalisation integral of 3-body decays */

//********************************************************************************
//Function protocols.

double Edndp3 (double yr, double ptr, double phirin, int res_num);
double norm3int (double x, void *paranorm);
double dnpir2N (double phi, void *para1);
double dnpir1N (double costh, void *para1);
double dn2ptN (double w2, void *para1);
double dn3ptN (double x, void *para1);
double Edndp3_2bodyN (double y, double pt, double phi, double m1, double m2,
		      double mr, int res_num);
double Edndp3_3bodyN (double y, double pt, double phi, double m1, double m2,
		      double m3, double mr, double norm3, int res_num);


/********************************************************************
*********************************************************************
*
*		decay routines
*
*********************************************************************
*********************************************************************/

double
norm3int (double x, void *paranorm)
{
  nblock *tmp = (nblock *) paranorm;
  double res = sqrt ((tmp->a - x) * (tmp->b - x)
		     * (x - tmp->c) * (x - tmp->d)) / x;
  return res;
}

//***********************************************************************
//The Main integrand for integration.

double dnpir2N (double phi, void *para1)
{
  pblockN *para = (pblockN *) para1;
  double D;
  double eR, plR, ptR, yR, phiR, sume, jac;
  double cphiR, sphiR;
  double dnr;			/* dn/mtdmt of resonance */


  sume = para->e + para->e0;

  D = para->e * para->e0 + para->pl * para->p0 * para->costh +
    para->pt * para->p0 * para->sinth * cos (phi) + para->m1 * para->m1;

  eR = para->mr * (sume * sume / D - 1.0);
  jac = para->mr + eR;
  plR = para->mr * sume * (para->pl - para->p0 * para->costh) / D;
  ptR = (eR * eR - plR * plR - para->mr * para->mr);

  if (ptR < 0.0)
    ptR = 0.0;

  else
    ptR = sqrt (ptR);

  yR = 0.5 * log ((eR + plR) / (eR - plR));
  cphiR = -jac * (para->p0 * para->sinth * cos (phi + para->phi)
		  - para->pt * cos (para->phi)) / (sume * ptR);
  sphiR = -jac * (para->p0 * para->sinth * sin (phi + para->phi)
		  - para->pt * sin (para->phi)) / (sume * ptR);

  if ((fabs (cphiR) > 1.000) || (fabs (sphiR) > 1.000))
    {
      if ((fabs (cphiR) > 1.0001) || (fabs (sphiR) > 1.0001))
	{
	  //printf ("  |phir| = %15.8lf  > 1 ! \n", phiR);
	  printf (" phi %15.8le D %15.8le \n", phi, D);
	  printf (" eR %15.8le plR %15.8le \n", eR, plR);
	  printf (" ptR %15.8le jac %15.8le \n", ptR, jac);
	  printf (" sume %15.8le costh %15.8le \n", sume, para->costh);

	  printf (" pt %15.8le \n", para->pt);
	  printf (" mt  %15.8le \n", para->mt);
	  printf (" y %15.8le \n", para->y);
	  printf (" e %15.8le \n", para->e);
	  printf (" e0 %15.8le \n", para->e0);
	  printf (" p0 %15.8le \n", para->p0);
	  printf (" pl %15.8le \n", para->pl);
	  printf (" phi %15.8le \n", para->phi);

	  printf (" m1 %15.8le \n", para->m1);
	  printf (" m2 %15.8le \n", para->m2);
	  printf (" m3 %15.8le \n", para->m3);
	  printf (" mr %15.8le \n", para->mr);
	  exit (0);
	}
      else
	{
	  if (cphiR > 1.0)
	    cphiR = 1.0;
	  if (cphiR < -1.0)
	    cphiR = -1.0;
	}
    }

  phiR = acos (cphiR);
  if (sphiR < 0.0)
    phiR = 2.0 * PI - phiR;
/*  if(isnan(yR)) printf("dNpir2N yR \n");
  if(isnan(ptR)) printf("dNpir2N ptR \n");
  if(isnan(phiR)) printf("dNpir2N phiR \n");
  if(isnan(para->res_num)) printf("dNpir2N res_num \n");*/
  dnr = Edndp3 (yR, ptR, phiR, para->res_num);

  /*printf(" phir = %15.8lf  ! ", phiR);
     printf(" ptR %15.8le jac %15.8le ", ptR, jac );
     printf(" dnr %15.8le \n", dnr); */
  return dnr * jac * jac / (2.0 * sume * sume);
}


double dnpir1N (double costh, void *para1)
{
  pblockN *para = (pblockN *) para1;
  double r;
  para->costh = costh;
  para->sinth = sqrt (1.0 - para->costh * para->costh);
  r = gauss (PTN2, *dnpir2N, 0.0, 2.0 * PI, para); //Integrates the "dnpir2N" kernal over phi using gaussian integration
  if(isnan(r))
  {
      printf("dnpir1N r \n");
  }
  return r;
}

double dn2ptN (double w2, void *para1)
{
  double test;
  pblockN *para = (pblockN *) para1;
  para->e0 = (para->mr * para->mr + para->m1 * para->m1 - w2) / (2 * para->mr);
  para->p0 = sqrt (para->e0 * para->e0 - para->m1 * para->m1);
/*  if(isnan(para->p0))
  {
      printf("dn2ptN, %E %E %E \n", para->p0, para->e0, para->m1);
  }*/
//   test = gauss (PTN1, *dnpir1N, -1.0, 1.0, para);
//   if(isnan(test))
//   {
//       printf("gauss, %E \n", test);
//   }
  return gauss (PTN1, *dnpir1N, -1.0, 1.0, para); //Integrate the "dnpir1N" kernal over cos(theta) using gaussian integration
}


double dn3ptN (double x, void* para1) //The integration kernal for "W" in 3-body decays.
{
  pblockN *para = (pblockN *) para1;
  double e0 =(para->mr * para->mr + para->m1 * para->m1 - x) / (2 * para->mr);
  double p0 = sqrt (e0 * e0 - para->m1 * para->m1);
  double a = (para->m2 + para->m3) * (para->m2 + para->m3);
  double b = (para->m2 - para->m3) * (para->m2 - para->m3);
  double re = p0 * sqrt ((x - a) * (x - b)) / x * dn2ptN (x, para);
  return re;
}


/********************************************************************
*
*	Edndp3_2bodyN()
*
* transverse momentum spectrum in GeV^-2 from pions out of resonances
*********************************************************************/
double Edndp3_2bodyN (double y, double pt, double phi, double m1, double m2,
		      double mr, int res_num)
		// /* in units of GeV^-2,includes phasespace and volume,
		//    does not include degeneracy factors  */
  //    double y;			/* rapidity of particle 1       */
  //    double pt;			/* transverse momentum of particle 1    */
  //    double phi;		/* phi angle of particle 1      */
  //    double m1, m2;		/* restmasses of decay particles in MeV */
  //    double mr;			/* restmass of resonance MeV            */
  //    int res_num;		/* Montecarlo number of the Resonance   */
{
  double mt = sqrt (pt * pt + m1 * m1);
  double norm2;			/* 2-body normalization         */
  pblockN para;
  double res2;

  para.pt = pt;
  para.mt = mt;
  para.e = mt * cosh (y);
  para.pl = mt * sinh (y);
  para.y = y;
  para.phi = phi;
  para.m1 = m1;
  para.m2 = m2;
  para.mr = mr;

  para.res_num = res_num;

  norm2 = 1.0 / (2.0 * PI);
  res2 = norm2 * dn2ptN (m2 * m2, &para); //Calls the integration routines for 2-body
/*  if(isnan(res2))
  {
      printf("res2, %E %E \n",res2, dn2ptN(m2*m2, &para));
  }*/
  return res2;			/* like Ed3ndp3_2body() */
}


/********************************************************************
*
*	Edndp3_3bodyN()
*
* transverse momentum spectrum in GeV^-2 from pions out of resonances
*********************************************************************/
double
Edndp3_3bodyN (double y, double pt, double phi, double m1, double m2,
	       double m3, double mr, double norm3, int res_num)
		/* in units of GeV^-2,includes phasespace and volume,
		   does not include degeneracy factors  */
{
  double mt = sqrt (pt * pt + m1 * m1);
  pblockN para;
  double wmin, wmax;
  double res3;
  double slope;			/* slope of resonance for high mt */
  int pn;

  para.pt = pt;
  para.mt = mt;
  para.y = y;
  para.e = mt * cosh (y);
  para.pl = mt * sinh (y);
  para.phi = phi;

  para.m1 = m1;
  para.m2 = m2;
  para.m3 = m3;
  para.mr = mr;

  pn = partid[MHALF + res_num];

  para.res_num = res_num;

  wmin = (m2 + m3) * (m2 + m3);
  wmax = (mr - m1) * (mr - m1);
  res3 = 2.0 * norm3 * gauss (PTS4, *dn3ptN, wmin, wmax, &para) / mr;  //Integrates "W" using gaussian
  return res3;
}



/**************************************************************************
*									  *
*	add_reso							  *
*									  *
* computes the pt, mt distribution including resonance decays		  *
***************************************************************************/


void add_reso (int pn, int pnR, int k, int j)
     // int pn;			/* internal number of daughter part. */
     // int pnR;			/* internal number of Resonance */
     // int k;			/* daughter-particle place in the decay
				 //   array  */
     // int j;			/* internal number of decay */
{
  nblock paranorm;		/* for 3body normalization integral */
  double y = 0.0;
  double m1, m2, m3, mr;
  double norm3;			/* normalisation of 3-body integral */
  int pn2, pn3, pn4;		/* internal numbers for resonances */
  int part;
  int l, i;
  int npt, nphi;

  npt = particle[pn].npt;
  nphi = particle[pn].nphi;

  // Determine the number of particles involved in the decay with the switch
  switch (abs (particleDecay[j].numpart))
    {
    case 1: //Only 1 particle, if it gets here, by accident, this prevents any integration for 1 particle chains
      break;

    case 2: // 2-body decay
      {
	if (k == 0)
	  pn2 = partid[MHALF + particleDecay[j].part[1]];

	else
	  pn2 = partid[MHALF + particleDecay[j].part[0]];

	//printf ("case 2:  i %3i j %3i k %3i \n", pn, j, k);
	m1 = particle[pn].mass;
	m2 = particle[pn2].mass;
	mr = particle[pnR].mass;

	while ((m1 + m2) > mr)
	  {
	    mr += 0.25 * particle[pnR].width;
	    m1 -= 0.5 * particle[pn].width;
	    m2 -= 0.5 * particle[pn2].width;
	  }

	for (l = 0; l < npt; l++)
	  {
	    for (i = 0; i < nphi; i++)
	      {
		// Call the 2-body decay integral and add its contribution to the daughter particle of interest
		particle[pn].dNdptdphi[l][i] += particleDecay[j].branch *
		  Edndp3_2bodyN (y, particle[pn].pt[l], PHI[i],
				 m1, m2, mr, particle[pnR].monval);
	      }
	  }
	break;
      }

    case 3: //3-body decay
      {
	if (k == 0)
	  {
	    pn2 = partid[MHALF + particleDecay[j].part[1]];
	    pn3 = partid[MHALF + particleDecay[j].part[2]];
	  }
	else
	  {
	    if (k == 1)
	      {
		pn2 = partid[MHALF + particleDecay[j].part[0]];
		pn3 = partid[MHALF + particleDecay[j].part[2]];
	      }
	    else
	      {
		pn2 = partid[MHALF + particleDecay[j].part[0]];
		pn3 = partid[MHALF + particleDecay[j].part[1]];
	      }
	  }

	m1 = particle[pn].mass;
	m2 = particle[pn2].mass;
	m3 = particle[pn3].mass;
	mr = particle[pnR].mass;
	paranorm.a = (mr + m1) * (mr + m1);
	paranorm.b = (mr - m1) * (mr - m1);
	paranorm.c = (m2 + m3) * (m2 + m3);
	paranorm.d = (m2 - m3) * (m2 - m3);
	norm3 = mr * mr / (2 * PI * gauss (PTS3, norm3int, paranorm.c,
					   paranorm.b, &paranorm));

	// printf("case 3:  i %3i j %3i k %3i \n",pn,j,k);

	for (i = 0; i < nphi; i++)
	  {
	    for (l = 0; l < npt; l++)
	      {
		// Call the 3-body decay integral and add its contribution to the daughter particle of interest
		particle[pn].dNdptdphi[l][i] += particleDecay[j].branch *
		  Edndp3_3bodyN (y, particle[pn].pt[l], PHI[i],
				 m1, m2, m3, mr, norm3, particle[pnR].monval);
	      }
	  }
	break;
      }

    case 4: //4-body decay (rare and low contribution)
      {
	if (k == 0)
	  {
	    pn2 = partid[MHALF + particleDecay[j].part[1]];
	    pn3 = partid[MHALF + particleDecay[j].part[2]];
	    pn4 = partid[MHALF + particleDecay[j].part[3]];
	  }
	else
	  {
	    if (k == 1)
	      {
		pn2 = partid[MHALF + particleDecay[j].part[0]];
		pn3 = partid[MHALF + particleDecay[j].part[2]];
		pn4 = partid[MHALF + particleDecay[j].part[3]];
	      }
	    else
	      {
		if (k == 2)
		  {
		    pn2 = partid[MHALF + particleDecay[j].part[0]];
		    pn3 = partid[MHALF + particleDecay[j].part[1]];
		    pn4 = partid[MHALF + particleDecay[j].part[3]];
		  }
		else
		  {
		    pn2 = partid[MHALF + particleDecay[j].part[0]];
		    pn3 = partid[MHALF + particleDecay[j].part[1]];
		    pn4 = partid[MHALF + particleDecay[j].part[2]];
		  }
	      }
	  }
	//approximate the 4-body with a 3-body decay with the 4th particle being the center of mass of 2 particles.
	m1 = particle[pn].mass;
	m2 = particle[pn2].mass;
	mr = particle[pnR].mass;
	m3 = 0.5 * (particle[pn3].mass + particle[pn4].mass + mr - m1 - m2);
	paranorm.a = (mr + m1) * (mr + m1);
	paranorm.b = (mr - m1) * (mr - m1);
	paranorm.c = (m2 + m3) * (m2 + m3);
	paranorm.d = (m2 - m3) * (m2 - m3);
	norm3 = mr * mr / (2 * PI * gauss (PTS3, norm3int, paranorm.c,
					   paranorm.b, &paranorm));
	// printf("case 3:  i %3i j %3i k %3i \n",pn,j,k);


	for (i = 0; i < nphi; i++)
	  {
	    for (l = 0; l < npt; l++)
	      {
		// the 4-body particleDecay approximated by the 3-body decay routine
		particle[pn].dNdptdphi[l][i] += particleDecay[j].branch *
		  Edndp3_3bodyN (y, particle[pn].pt[l], PHI[i],
				 m1, m2, m3, mr, norm3, particle[pnR].monval);
	      }
	  }
	break;
      }

    default:
      printf ("ERROR in add_reso! \n");
      printf ("%i decay not implemented ! \n", abs (particleDecay[j].numpart));
      exit (0);
    }
}




/**************************************************************************
*									  *
*	cal_pt_spec							  *
*									  *
* conductr p_1440()s add_reso for computing the pt, mt distribution including	  *
* resonance decays		 					  *
***************************************************************************/


void cal_reso_decays (int maxpart, int maxdecay, int bound)
     // int maxpart, maxdecay;
     // int bound;			/* monte-carlo value of a lower bound
				 //   in the spectrum of particles taken
				 //   into  account                      */
{
  int i, j, k, l, ll;
  int pn, pnR, pnaR;
  int part;

  printf (" CALCULATE RESONANCE DECAYS (as fast as I can) \n");
  pn = partid[MHALF + bound];
  for(i=maxpart-1;i > pn-1;i--)  //Cycle the particles known from the resoweak.dat input
    {
      part = particle[i].monval;
      printf ("Calculating the decays with ");
      printf ("%s \n", particle[i].name);
      printf ("%i, %i", part, maxdecay);
      switch (particle[i].baryon)  // Check to see whether or not the particle is baryon, anti-baryon or meson
	{
	case 1: //Baryon
	  {
	    printf("Is a baryon. \n");

	    for (j = 0; j < maxdecay; j++) // Cycle through every decay channel known (as given in resoweak.dat)
	      {                            // to see if the particle was a daughter particle in a decay channel
		pnR = partid[MHALF + particleDecay[j].reso];
		//printf("Partid is %i.\n",pnR);
		for (k = 0; k < abs (particleDecay[j].numpart); k++)
		  {
		    if ((part == particleDecay[j].part[k]) && (particleDecay[j].numpart != 1))// Make sure that the particleDecay channel isn't trivial
		      {                                                       // and contains the daughter particle
			//printf("Calculating a decay \n");
			add_reso (i, pnR, k, j);
		      }
		  }
	      }
	    break;
	  }

	case -1: //Anti-Baryon
	  {
	    printf("Is an anti-baryon.\n");
	    for (j = 0; j < maxdecay; j++)// Cycle through every decay channel known (as given in resoweak.dat)
	      {                            // to see if the particle was a daughter particle in a decay channel
		pnaR = partid[MHALF - particleDecay[j].reso];
		//printf("Partid is %i.\n",pnaR);
		for (k = 0; k < abs (particleDecay[j].numpart); k++)
		  {
		    if ((-part == particleDecay[j].part[k]) && (particleDecay[j].numpart != 1))// Make sure that the decay channel isn't trivial
		      {                                                        // and contains the daughter particle
			//printf("Calculating a decay \n");
			add_reso (i, pnaR, k, j);
		      }
		  }
	      }
	    break;
	  }

	case 0:// Meson
	  {
	    printf("Is a meson. \n");

	    for (j = 0; j < maxdecay; j++)
	      {
		pnR = partid[MHALF + particleDecay[j].reso];
		//printf("Partid is %i.\n",pnR);
		for (k = 0; k < abs (particleDecay[j].numpart); k++)
		  {
		    if (particle[pnR].baryon == 1)
		      {
			pnaR = partid[MHALF - particleDecay[j].reso];
			if ((particle[i].charge == 0)
			    && (particle[i].strange == 0))
			  {
			    if ((part == particleDecay[j].part[k])
				&& (particleDecay[j].numpart != 1))
			      {
				//printf("Calculating a decay \n");
				add_reso (i, pnR, k, j);
				add_reso (i, pnaR, k, j);
			      }
			  }
			else
			  {
			    if ((part == particleDecay[j].part[k])
				&& (particleDecay[j].numpart != 1))
			      {
				//printf("Calculating a decay \n");
				add_reso (i, pnR, k, j);
			      }
			    if ((-part == particleDecay[j].part[k])
				&& (particleDecay[j].numpart != 1))
			      {
				//printf("Calculating a decay \n");
				add_reso (i, pnaR, k, j);
			      }
			  }
		      }
		    else
		      {
			if ((part == particleDecay[j].part[k])
			    && (particleDecay[j].numpart != 1))
			  {
			    //printf("Calculating a decay \n");
			    add_reso (i, pnR, k, j);
			  }
		      }
		  }
	      }
	    break;
	  }


	default:
	  printf ("Error in switch in func partden_wdecay \n");
	  exit (0);

	}
    }

  //***********************************************************************************************
  // Routine to look at individual resonance contributions.  Must set the above loop to END right before
  // the particle.  Example, i=3 is pi(-).  Must stop the above loop to stop before i=3 to examine decays
  // containing the pi(-) as a daughter

  /*

  int test = 0;
  int testmax = 9; // The # of resonances to look at.
  while(test < testmax)
    {
      i = 3;
      //for (i = maxpart - 1; i >= maxpart -1; i--)
      part = particle[i].monval;
      printf ("Calculating the decays for ");
      printf ("%s \n", particle[i].name);
      //printf ("%i, %i", part, maxdecay);
      switch (particle[i].baryon)
	{
	case 1:
	  {
	    //printf("Is a baryon. \n");

	    for (j = 0; j < maxdecay; j++)
	      {
		pnR = partid[MHALF + decay[j].reso];
		//printf("Partid is %i.\n",pnR);
		for (k = 0; k < abs (decay[j].numpart); k++)
		  {
		    if ((part == decay[j].part[k]) && (decay[j].numpart != 1))
		      {
			if (decay[j].numpart == 2)
			  printf("Calculating decay %i -> %i + %i \n", decay[j].reso, decay[j].part[0], decay[j].part[1]);
			if (decay[j].numpart == 3)
			  printf("Calculating decay %i -> %i + %i + %i\n", decay[j].reso, decay[j].part[0], decay[j].part[1], decay[j].part[2]);
			if (decay[j].numpart == 4)
			  printf("Calculating decay %i -> %i + %i + %i + %i \n", decay[j].reso, decay[j].part[0], decay[j].part[1], decay[j].part[2], decay[j].part[3]);
			test++;
			add_reso (i, pnR, k, j);
		      }
		  }
	      }
	    break;
	  }

	case -1:
	  {
	    //printf("Is an anti-baryon.\n");
	    for (j = 0; j < maxdecay; j++)
	      {
		pnaR = partid[MHALF - decay[j].reso];
		//printf("Partid is %i.\n",pnaR);
		for (k = 0; k < abs (decay[j].numpart); k++)
		  {
		    if ((-part == decay[j].part[k]) && (decay[j].numpart != 1))
		      {
			if (decay[j].numpart == 2)
			  printf("Calculating decay %i -> %i + %i \n", decay[j].reso, decay[j].part[0], decay[j].part[1]);
			if (decay[j].numpart == 3)
			  printf("Calculating decay %i -> %i + %i + %i\n", decay[j].reso, decay[j].part[0], decay[j].part[1], decay[j].part[2]);
			if (decay[j].numpart == 4)
			  printf("Calculating decay %i -> %i + %i + %i + %i \n", decay[j].reso, decay[j].part[0], decay[j].part[1], decay[j].part[2], decay[j].part[3]);
			test++;
			add_reso (i, pnaR, k, j);
		      }
		  }
	      }
	    break;
	  }

	case 0:
	  {
	    //printf("Is a meson. \n");

	    for (j = 0; j < maxdecay; j++)
	      {
		pnR = partid[MHALF + decay[j].reso];
		//printf("Partid is %i.\n",pnR);
		for (k = 0; k < abs (decay[j].numpart); k++)
		  {
		    if (particle[pnR].baryon == 1 && test < testmax)
		      {
			pnaR = partid[MHALF - decay[j].reso];
			if ((particle[i].charge == 0)
			    && (particle[i].strange == 0))
			  {
			    if ((part == decay[j].part[k])
				&& (decay[j].numpart != 1 && test < testmax))
			      {

				if (decay[j].numpart == 2)
				  printf("Calculating decay %i -> %i + %i \n", decay[j].reso, decay[j].part[0], decay[j].part[1]);
				if (decay[j].numpart == 3)
				  printf("Calculating decay %i -> %i + %i + %i\n", decay[j].reso, decay[j].part[0], decay[j].part[1], decay[j].part[2]);
				if (decay[j].numpart == 4)
				  printf("Calculating decay %i -> %i + %i + %i + %i \n", decay[j].reso, decay[j].part[0], decay[j].part[1], decay[j].part[2], decay[j].part[3]);
				test++;
				add_reso (i, pnR, k, j);
				add_reso (i, pnaR, k, j);
			      }
			  }
			else
			  {
			    if ((part == decay[j].part[k])
				&& (decay[j].numpart != 1) && test < testmax)
			      {
				if (decay[j].numpart == 2)
				  printf("Calculating decay %i -> %i + %i \n", decay[j].reso, decay[j].part[0], decay[j].part[1]);
				if (decay[j].numpart == 3)
				  printf("Calculating decay %i -> %i + %i + %i\n", decay[j].reso, decay[j].part[0], decay[j].part[1], decay[j].part[2]);
				if (decay[j].numpart == 4)
				  printf("Calculating decay %i -> %i + %i + %i + %i \n", decay[j].reso, decay[j].part[0], decay[j].part[1], decay[j].part[2], decay[j].part[3]);
				test++;
				add_reso (i, pnR, k, j);
			      }
			    if ((-part == decay[j].part[k])
				&& (decay[j].numpart != 1 && test < testmax))
			      {
				if (decay[j].numpart == 2)
				  printf("Calculating decay %i -> %i + %i \n", decay[j].reso, decay[j].part[0], decay[j].part[1]);
				if (decay[j].numpart == 3)
				  printf("Calculating decay %i -> %i + %i + %i\n", decay[j].reso, decay[j].part[0], decay[j].part[1], decay[j].part[2]);
				if (decay[j].numpart == 4)
				  printf("Calculating decay %i -> %i + %i + %i + %i \n", decay[j].reso, decay[j].part[0], decay[j].part[1], decay[j].part[2], decay[j].part[3]);
				test++;
				add_reso (i, pnaR, k, j);
			      }
			  }
		      }
		    else
		      {
			if ((part == decay[j].part[k])
			    && (decay[j].numpart != 1 && test < testmax))
			  {
			    if (decay[j].numpart == 2)
			      printf("Calculating decay %i -> %i + %i \n", decay[j].reso, decay[j].part[0], decay[j].part[1]);
			    if (decay[j].numpart == 3)
			      printf("Calculating decay %i -> %i + %i + %i\n", decay[j].reso, decay[j].part[0], decay[j].part[1], decay[j].part[2]);
			    if (decay[j].numpart == 4)
			      printf("Calculating decay %i -> %i + %i + %i + %i \n", decay[j].reso, decay[j].part[0], decay[j].part[1], decay[j].part[2], decay[j].part[3]);
			    test++;
			    add_reso (i, pnR, k, j);
			  }
		      }
		  }
	      }
	    break;
	  }


	default:
	  printf ("Error in switch in func partden_wdecay \n");
	  exit (0);

	}
    }

  */
}
