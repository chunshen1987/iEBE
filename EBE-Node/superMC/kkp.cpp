/* kkp.f -- translated by f2c (version 20090411).
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/


/* ===================================================================== */

/*     ------------------------------------------------------------ */
/*     Fragmentation functions for: Pions, Kaons, Protons, Neutrons */
/*            (includes mass-threshholds for c and b quarks) */

/*     Reference: B.A.Kniehl, G.Kramer, B.Potter, NPB582 (2000) 514 */
/*     ------------------------------------------------------------ */

/*     ih, iset, x, qs are input; dh is output. */
/*     ih   = 1 : (pi^+ + pi^-)  /2 */
/*     ih   = 2 : (K^+ + K^-)    /2 */
/*     ih   = 3 : (K^0 + K^0_bar)/2 */
/*     ih   = 4 : (p + p_bar)    /2 */
/*     ih   = 5 : (pi^0) */
/*     ih   = 6 : (n + n_bar)    /2 */
/*     ih   = 7 : (h^+ + h^-)         [as sum of pions, kaons and protons] */

/*     iset = 0 : LO */
/*     iset = 1 : NLO */

/*     x    = longitudinal-momentum fraction */
/*     qs   = fragmentation scale (in GeV) */

/*     Parton label: */
/*     0    1    2    3    4    5    6    7    8     9    10 */
/*     g    u   ubar  d   dbar  s   sbar  c   cbar   b   bbar */

/*     Lambda_QCD (in GeV): */
/*     0.088 in LO */
/*     0.213 in NLO */

/* ===================================================================== */
void KLNModel::kkpFF(int *ih, int *iset, double *x, double *qs, double *dh)
{
    /* System generated locals */
    double d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9, d__10, 
	    d__11;

    /* Local variables */
    static double a[11], b[3], s, q0, sb, sc, dkb, dkc, dkd, dpb, dpc, 
	    dkg, dpd, dpg, dks, dku, dps, dpu, dk0b, dk0c, dk0d, dk0s, rmbb, 
	    rmcc, dprb, dprc, dprd, rlam, dprg, dprs, dpru;

/* --- Mass-thresholds: */
    rmcc = 2.9788;
    rmbb = 9.46037;
/* --- Q_0 (in GeV): */
    q0 = sqrt(2.);
/* --- BAK */
    if (*qs < q0) {
	*qs = q0;
    }
/* --- LO FFs */
    if (*iset == 0) {
	rlam = .088;
/* Computing 2nd power */
	d__1 = *qs;
/* Computing 2nd power */
	d__2 = rlam;
/* Computing 2nd power */
	d__3 = q0;
/* Computing 2nd power */
	d__4 = rlam;
	s = log(log(d__1 * d__1 / (d__2 * d__2)) / log(d__3 * d__3 / (d__4 * 
		d__4)));
/* Computing 2nd power */
	d__1 = *qs;
/* Computing 2nd power */
	d__2 = rlam;
/* Computing 2nd power */
	d__3 = rmcc;
/* Computing 2nd power */
	d__4 = rlam;
	sc = log(log(d__1 * d__1 / (d__2 * d__2)) / log(d__3 * d__3 / (d__4 * 
		d__4)));
/* Computing 2nd power */
	d__1 = *qs;
/* Computing 2nd power */
	d__2 = rlam;
/* Computing 2nd power */
	d__3 = rmbb;
/* Computing 2nd power */
	d__4 = rlam;
	sb = log(log(d__1 * d__1 / (d__2 * d__2)) / log(d__3 * d__3 / (d__4 * 
		d__4)));
/* ---------------------- LO PION ------------------------------ */
	b[0] = 6.0451;
	b[1] = -.71378;
	b[2] = 2.92133;
	a[0] = -6.61523;
	a[1] = -1.64978;
	a[2] = 2.68223;
	a[3] = .14705;
	a[4] = -1.08423;
	a[5] = -.43182;
	a[6] = 1.48429;
	a[7] = 1.32887;
	a[8] = -1.78696;
	a[9] = .23086;
	a[10] = -.29182;
/* Computing 2nd power */
	d__1 = s;
/* Computing 3rd power */
	d__2 = s;
/* Computing 2nd power */
	d__4 = s;
/* Computing 3rd power */
	d__5 = s;
	d__3 = b[1] + a[3] * s + a[4] * (d__4 * d__4) + a[5] * (d__5 * (d__5 *
		 d__5));
	d__6 = 1. - *x;
/* Computing 2nd power */
	d__8 = s;
/* Computing 3rd power */
	d__9 = s;
	d__7 = b[2] + a[6] * s + a[7] * (d__8 * d__8) + a[8] * (d__9 * (d__9 *
		 d__9));
/* Computing 2nd power */
	d__10 = s;
	dpg = (b[0] + a[0] * s + a[1] * (d__1 * d__1) + a[2] * (d__2 * (d__2 *
		 d__2))) * pow__dd(x, &d__3) * pow__dd(&d__6, &d__7) * ((a[9] * 
		s + a[10] * (d__10 * d__10)) / *x + 1.);
	b[0] = .5461;
	b[1] = -1.46616;
	b[2] = 1.01864;
	a[0] = -.22946;
	a[1] = -.22594;
	a[2] = .21119;
	a[3] = -.45404;
	a[4] = -.12684;
	a[5] = .27646;
	a[6] = .95367;
	a[7] = -1.09835;
	a[8] = .74657;
	a[9] = -.01877;
	a[10] = .02949;
/* Computing 2nd power */
	d__1 = s;
/* Computing 3rd power */
	d__2 = s;
/* Computing 2nd power */
	d__4 = s;
/* Computing 3rd power */
	d__5 = s;
	d__3 = b[1] + a[3] * s + a[4] * (d__4 * d__4) + a[5] * (d__5 * (d__5 *
		 d__5));
	d__6 = 1. - *x;
/* Computing 2nd power */
	d__8 = s;
/* Computing 3rd power */
	d__9 = s;
	d__7 = b[2] + a[6] * s + a[7] * (d__8 * d__8) + a[8] * (d__9 * (d__9 *
		 d__9));
/* Computing 2nd power */
	d__10 = s;
	dpu = (b[0] + a[0] * s + a[1] * (d__1 * d__1) + a[2] * (d__2 * (d__2 *
		 d__2))) * pow__dd(x, &d__3) * pow__dd(&d__6, &d__7) * ((a[9] * 
		s + a[10] * (d__10 * d__10)) / *x + 1.);
	b[0] = 22.2815;
	b[1] = .12732;
	b[2] = 6.13697;
	a[0] = -20.8125;
	a[1] = -11.5725;
	a[2] = 15.5372;
	a[3] = .23075;
	a[4] = -2.71424;
	a[5] = 1.72456;
	a[6] = 2.18849;
	a[7] = -5.04475;
	a[8] = 3.29117;
	a[9] = .09044;
	a[10] = -.07589;
/* Computing 2nd power */
	d__1 = s;
/* Computing 3rd power */
	d__2 = s;
/* Computing 2nd power */
	d__4 = s;
/* Computing 3rd power */
	d__5 = s;
	d__3 = b[1] + a[3] * s + a[4] * (d__4 * d__4) + a[5] * (d__5 * (d__5 *
		 d__5));
	d__6 = 1. - *x;
/* Computing 2nd power */
	d__8 = s;
/* Computing 3rd power */
	d__9 = s;
	d__7 = b[2] + a[6] * s + a[7] * (d__8 * d__8) + a[8] * (d__9 * (d__9 *
		 d__9));
/* Computing 2nd power */
	d__10 = s;
	dps = (b[0] + a[0] * s + a[1] * (d__1 * d__1) + a[2] * (d__2 * (d__2 *
		 d__2))) * pow__dd(x, &d__3) * pow__dd(&d__6, &d__7) * ((a[9] * 
		s + a[10] * (d__10 * d__10)) / *x + 1.);
	b[0] = 8.755;
	b[1] = -.38611;
	b[2] = 5.61846;
	a[0] = -9.32277;
	a[1] = 1.806;
	a[2] = 2.02179;
	a[3] = -.4119;
	a[4] = -.48496;
	a[5] = .42525;
	a[6] = .74035;
	a[7] = -.64929;
	a[8] = .66788;
	a[9] = .06652;
	a[10] = -.05531;
/* Computing 2nd power */
	d__1 = sc;
/* Computing 3rd power */
	d__2 = sc;
/* Computing 2nd power */
	d__4 = sc;
/* Computing 3rd power */
	d__5 = sc;
	d__3 = b[1] + a[3] * sc + a[4] * (d__4 * d__4) + a[5] * (d__5 * (d__5 
		* d__5));
	d__6 = 1. - *x;
/* Computing 2nd power */
	d__8 = sc;
/* Computing 3rd power */
	d__9 = sc;
	d__7 = b[2] + a[6] * sc + a[7] * (d__8 * d__8) + a[8] * (d__9 * (d__9 
		* d__9));
/* Computing 2nd power */
	d__10 = sc;
	dpc = (b[0] + a[0] * sc + a[1] * (d__1 * d__1) + a[2] * (d__2 * (d__2 
		* d__2))) * pow__dd(x, &d__3) * pow__dd(&d__6, &d__7) * ((a[9] *
		 sc + a[10] * (d__10 * d__10)) / *x + 1.);
	b[0] = .31147;
	b[1] = -1.92993;
	b[2] = 3.47086;
	a[0] = -.19319;
	a[1] = -.10487;
	a[2] = .18824;
	a[3] = -.44692;
	a[4] = -.08271;
	a[5] = .30441;
	a[6] = .79775;
	a[7] = -.28091;
	a[8] = .39504;
	a[9] = -.04887;
	a[10] = .03212;
/* Computing 2nd power */
	d__1 = sb;
/* Computing 3rd power */
	d__2 = sb;
/* Computing 2nd power */
	d__4 = sb;
/* Computing 3rd power */
	d__5 = sb;
	d__3 = b[1] + a[3] * sb + a[4] * (d__4 * d__4) + a[5] * (d__5 * (d__5 
		* d__5));
	d__6 = 1. - *x;
/* Computing 2nd power */
	d__8 = sb;
/* Computing 3rd power */
	d__9 = sb;
	d__7 = b[2] + a[6] * sb + a[7] * (d__8 * d__8) + a[8] * (d__9 * (d__9 
		* d__9));
/* Computing 2nd power */
	d__10 = sb;
	dpb = (b[0] + a[0] * sb + a[1] * (d__1 * d__1) + a[2] * (d__2 * (d__2 
		* d__2))) * pow__dd(x, &d__3) * pow__dd(&d__6, &d__7) * ((a[9] *
		 sb + a[10] * (d__10 * d__10)) / *x + 1.);
/* ---------------------- LO KAON ------------------------------ */
	b[0] = .02862;
	b[1] = -2.94091;
	b[2] = 2.73474;
	a[0] = -.02113;
	a[1] = .00389;
	a[2] = .00901;
	a[3] = .66881;
	a[4] = -.2967;
	a[5] = .20574;
	a[6] = -.58222;
	a[7] = .04329;
	a[8] = .78033;
	a[9] = .03586;
	a[10] = -.0122;
/* Computing 2nd power */
	d__1 = s;
/* Computing 3rd power */
	d__2 = s;
/* Computing 2nd power */
	d__4 = s;
/* Computing 3rd power */
	d__5 = s;
	d__3 = b[1] + a[3] * s + a[4] * (d__4 * d__4) + a[5] * (d__5 * (d__5 *
		 d__5));
	d__6 = 1. - *x;
/* Computing 2nd power */
	d__8 = s;
/* Computing 3rd power */
	d__9 = s;
	d__7 = b[2] + a[6] * s + a[7] * (d__8 * d__8) + a[8] * (d__9 * (d__9 *
		 d__9));
/* Computing 2nd power */
	d__10 = s;
	dkg = (b[0] + a[0] * s + a[1] * (d__1 * d__1) + a[2] * (d__2 * (d__2 *
		 d__2))) * pow__dd(x, &d__3) * pow__dd(&d__6, &d__7) * ((a[9] * 
		s + a[10] * (d__10 * d__10)) / *x + 1.);
	b[0] = .25937;
	b[1] = -.61925;
	b[2] = .85946;
	a[0] = -.10502;
	a[1] = .00572;
	a[2] = -.00269;
	a[3] = .09956;
	a[4] = .07389;
	a[5] = -7e-4;
	a[6] = .57965;
	a[7] = .26397;
	a[8] = -.12764;
	a[9] = .15303;
	a[10] = .14807;
/* Computing 2nd power */
	d__1 = s;
/* Computing 3rd power */
	d__2 = s;
/* Computing 2nd power */
	d__4 = s;
/* Computing 3rd power */
	d__5 = s;
	d__3 = b[1] + a[3] * s + a[4] * (d__4 * d__4) + a[5] * (d__5 * (d__5 *
		 d__5));
	d__6 = 1. - *x;
/* Computing 2nd power */
	d__8 = s;
/* Computing 3rd power */
	d__9 = s;
	d__7 = b[2] + a[6] * s + a[7] * (d__8 * d__8) + a[8] * (d__9 * (d__9 *
		 d__9));
/* Computing 2nd power */
	d__10 = s;
	dku = (b[0] + a[0] * s + a[1] * (d__1 * d__1) + a[2] * (d__2 * (d__2 *
		 d__2))) * pow__dd(x, &d__3) * pow__dd(&d__6, &d__7) * ((a[9] * 
		s + a[10] * (d__10 * d__10)) / *x + 1.);
	b[0] = 5.38115;
	b[1] = -.00321;
	b[2] = 3.07632;
	a[0] = -3.05084;
	a[1] = -1.10056;
	a[2] = 1.31207;
	a[3] = -.25889;
	a[4] = -.18494;
	a[5] = .13994;
	a[6] = 1.13745;
	a[7] = -.90413;
	a[8] = .56581;
	a[9] = .05141;
	a[10] = -.00697;
/* Computing 2nd power */
	d__1 = s;
/* Computing 3rd power */
	d__2 = s;
/* Computing 2nd power */
	d__4 = s;
/* Computing 3rd power */
	d__5 = s;
	d__3 = b[1] + a[3] * s + a[4] * (d__4 * d__4) + a[5] * (d__5 * (d__5 *
		 d__5));
	d__6 = 1. - *x;
/* Computing 2nd power */
	d__8 = s;
/* Computing 3rd power */
	d__9 = s;
	d__7 = b[2] + a[6] * s + a[7] * (d__8 * d__8) + a[8] * (d__9 * (d__9 *
		 d__9));
/* Computing 2nd power */
	d__10 = s;
	dkd = (b[0] + a[0] * s + a[1] * (d__1 * d__1) + a[2] * (d__2 * (d__2 *
		 d__2))) * pow__dd(x, &d__3) * pow__dd(&d__6, &d__7) * ((a[9] * 
		s + a[10] * (d__10 * d__10)) / *x + 1.);
	b[0] = 5.18266;
	b[1] = -.17751;
	b[2] = 4.30306;
	a[0] = -3.48519;
	a[1] = -1.00982;
	a[2] = 1.17996;
	a[3] = .02309;
	a[4] = -.61327;
	a[5] = -.03532;
	a[6] = 1.00547;
	a[7] = -.51779;
	a[8] = .20683;
	a[9] = .13514;
	a[10] = -.17778;
/* Computing 2nd power */
	d__1 = sc;
/* Computing 3rd power */
	d__2 = sc;
/* Computing 2nd power */
	d__4 = sc;
/* Computing 3rd power */
	d__5 = sc;
	d__3 = b[1] + a[3] * sc + a[4] * (d__4 * d__4) + a[5] * (d__5 * (d__5 
		* d__5));
	d__6 = 1. - *x;
/* Computing 2nd power */
	d__8 = sc;
/* Computing 3rd power */
	d__9 = sc;
	d__7 = b[2] + a[6] * sc + a[7] * (d__8 * d__8) + a[8] * (d__9 * (d__9 
		* d__9));
/* Computing 2nd power */
	d__10 = sc;
	dkc = (b[0] + a[0] * sc + a[1] * (d__1 * d__1) + a[2] * (d__2 * (d__2 
		* d__2))) * pow__dd(x, &d__3) * pow__dd(&d__6, &d__7) * ((a[9] *
		 sc + a[10] * (d__10 * d__10)) / *x + 1.);
	b[0] = 1.57044;
	b[1] = -.84143;
	b[2] = 6.01488;
	a[0] = -1.7834;
	a[1] = .571;
	a[2] = .15469;
	a[3] = -.43448;
	a[4] = -.05314;
	a[5] = -.36621;
	a[6] = .72953;
	a[7] = -.64433;
	a[8] = .92351;
	a[9] = .01024;
	a[10] = -.0616;
/* Computing 2nd power */
	d__1 = sb;
/* Computing 3rd power */
	d__2 = sb;
/* Computing 2nd power */
	d__4 = sb;
/* Computing 3rd power */
	d__5 = sb;
	d__3 = b[1] + a[3] * sb + a[4] * (d__4 * d__4) + a[5] * (d__5 * (d__5 
		* d__5));
	d__6 = 1. - *x;
/* Computing 2nd power */
	d__8 = sb;
/* Computing 3rd power */
	d__9 = sb;
	d__7 = b[2] + a[6] * sb + a[7] * (d__8 * d__8) + a[8] * (d__9 * (d__9 
		* d__9));
/* Computing 2nd power */
	d__10 = sb;
	dkb = (b[0] + a[0] * sb + a[1] * (d__1 * d__1) + a[2] * (d__2 * (d__2 
		* d__2))) * pow__dd(x, &d__3) * pow__dd(&d__6, &d__7) * ((a[9] *
		 sb + a[10] * (d__10 * d__10)) / *x + 1.);
/* ---------------------- LO PROTON ----------------------------- */
	b[0] = .73953;
	b[1] = -.76986;
	b[2] = 7.69079;
	a[0] = -1.64519;
	a[1] = 1.01189;
	a[2] = -.10175;
	a[3] = -3.58787;
	a[4] = 13.8025;
	a[5] = -13.8902;
	a[6] = -2.8447;
	a[7] = -.36719;
	a[8] = -2.21825;
	a[9] = 1.26515;
	a[10] = -1.96117;
/* Computing 2nd power */
	d__1 = s;
/* Computing 3rd power */
	d__2 = s;
/* Computing 2nd power */
	d__4 = s;
/* Computing 3rd power */
	d__5 = s;
	d__3 = b[1] + a[3] * s + a[4] * (d__4 * d__4) + a[5] * (d__5 * (d__5 *
		 d__5));
	d__6 = 1. - *x;
/* Computing 2nd power */
	d__8 = s;
/* Computing 3rd power */
	d__9 = s;
	d__7 = b[2] + a[6] * s + a[7] * (d__8 * d__8) + a[8] * (d__9 * (d__9 *
		 d__9));
/* Computing 2nd power */
	d__10 = s;
/* Computing 3rd power */
	d__11 = s;
	dprg = (b[0] + a[0] * s + a[1] * (d__1 * d__1) + a[2] * (d__2 * (d__2 
		* d__2))) * pow__dd(x, &d__3) * pow__dd(&d__6, &d__7) * ((a[9] *
		 s + a[10] * (d__10 * d__10) + d__11 * (d__11 * d__11) * 
		.54769f) / *x + 1.);
	b[0] = .40211;
	b[1] = -.85973;
	b[2] = 2.8016;
	a[0] = -.21633;
	a[1] = -.07045;
	a[2] = .07831;
	a[3] = .13987;
	a[4] = -.82412;
	a[5] = .43114;
	a[6] = .78923;
	a[7] = -.05344;
	a[8] = .0146;
	a[9] = .05198;
	a[10] = -.04623;
/* Computing 2nd power */
	d__1 = s;
/* Computing 3rd power */
	d__2 = s;
/* Computing 2nd power */
	d__4 = s;
/* Computing 3rd power */
	d__5 = s;
	d__3 = b[1] + a[3] * s + a[4] * (d__4 * d__4) + a[5] * (d__5 * (d__5 *
		 d__5));
	d__6 = 1. - *x;
/* Computing 2nd power */
	d__8 = s;
/* Computing 3rd power */
	d__9 = s;
	d__7 = b[2] + a[6] * s + a[7] * (d__8 * d__8) + a[8] * (d__9 * (d__9 *
		 d__9));
/* Computing 2nd power */
	d__10 = s;
	dpru = (b[0] + a[0] * s + a[1] * (d__1 * d__1) + a[2] * (d__2 * (d__2 
		* d__2))) * pow__dd(x, &d__3) * pow__dd(&d__6, &d__7) * ((a[9] *
		 s + a[10] * (d__10 * d__10)) / *x + 1.);
	b[0] = 4.07885;
	b[1] = -.09735;
	b[2] = 4.99191;
	a[0] = -2.97392;
	a[1] = -.92973;
	a[2] = 1.23517;
	a[3] = .25834;
	a[4] = -1.52246;
	a[5] = .7706;
	a[6] = 1.14379;
	a[7] = -.8532;
	a[8] = .45607;
	a[9] = .07174;
	a[10] = -.08321;
/* Computing 2nd power */
	d__1 = s;
/* Computing 3rd power */
	d__2 = s;
/* Computing 2nd power */
	d__4 = s;
/* Computing 3rd power */
	d__5 = s;
	d__3 = b[1] + a[3] * s + a[4] * (d__4 * d__4) + a[5] * (d__5 * (d__5 *
		 d__5));
	d__6 = 1. - *x;
/* Computing 2nd power */
	d__8 = s;
/* Computing 3rd power */
	d__9 = s;
	d__7 = b[2] + a[6] * s + a[7] * (d__8 * d__8) + a[8] * (d__9 * (d__9 *
		 d__9));
/* Computing 2nd power */
	d__10 = s;
	dprs = (b[0] + a[0] * s + a[1] * (d__1 * d__1) + a[2] * (d__2 * (d__2 
		* d__2))) * pow__dd(x, &d__3) * pow__dd(&d__6, &d__7) * ((a[9] *
		 s + a[10] * (d__10 * d__10)) / *x + 1.);
	b[0] = .11061;
	b[1] = -1.5434;
	b[2] = 2.20681;
	a[0] = -.07726;
	a[1] = .05422;
	a[2] = -.03364;
	a[3] = -.20804;
	a[4] = .29038;
	a[5] = -.23662;
	a[6] = .62274;
	a[7] = .29713;
	a[8] = -.21861;
	a[9] = .00831;
	a[10] = 6.5e-4;
/* Computing 2nd power */
	d__1 = sc;
/* Computing 3rd power */
	d__2 = sc;
/* Computing 2nd power */
	d__4 = sc;
/* Computing 3rd power */
	d__5 = sc;
	d__3 = b[1] + a[3] * sc + a[4] * (d__4 * d__4) + a[5] * (d__5 * (d__5 
		* d__5));
	d__6 = 1. - *x;
/* Computing 2nd power */
	d__8 = sc;
/* Computing 3rd power */
	d__9 = sc;
	d__7 = b[2] + a[6] * sc + a[7] * (d__8 * d__8) + a[8] * (d__9 * (d__9 
		* d__9));
/* Computing 2nd power */
	d__10 = sc;
	dprc = (b[0] + a[0] * sc + a[1] * (d__1 * d__1) + a[2] * (d__2 * (
		d__2 * d__2))) * pow__dd(x, &d__3) * pow__dd(&d__6, &d__7) * ((
		a[9] * sc + a[10] * (d__10 * d__10)) / *x + 1.);
	b[0] = 40.0971;
	b[1] = .74249;
	b[2] = 12.3729;
	a[0] = -123.531;
	a[1] = 128.666;
	a[2] = -29.1808;
	a[3] = -1.29639;
	a[4] = -3.65003;
	a[5] = 3.0534;
	a[6] = -1.04932;
	a[7] = .34662;
	a[8] = -1.34412;
	a[9] = -.0429;
	a[10] = -.30359;
/* Computing 2nd power */
	d__1 = sb;
/* Computing 3rd power */
	d__2 = sb;
/* Computing 2nd power */
	d__4 = sb;
/* Computing 3rd power */
	d__5 = sb;
	d__3 = b[1] + a[3] * sb + a[4] * (d__4 * d__4) + a[5] * (d__5 * (d__5 
		* d__5));
	d__6 = 1. - *x;
/* Computing 2nd power */
	d__8 = sb;
/* Computing 3rd power */
	d__9 = sb;
	d__7 = b[2] + a[6] * sb + a[7] * (d__8 * d__8) + a[8] * (d__9 * (d__9 
		* d__9));
/* Computing 2nd power */
	d__10 = sb;
	dprb = (b[0] + a[0] * sb + a[1] * (d__1 * d__1) + a[2] * (d__2 * (
		d__2 * d__2))) * pow__dd(x, &d__3) * pow__dd(&d__6, &d__7) * ((
		a[9] * sb + a[10] * (d__10 * d__10)) / *x + 1.);
/* --- NLO FFs */
    } else {
	if (*iset != 1) {
	  cout << "ERRORin kkpFF: iset must be 0 (LO) or 1 (NLO)" << endl;
	  exit(1);
	}
	rlam = .213;
/* Computing 2nd power */
	d__1 = *qs;
/* Computing 2nd power */
	d__2 = rlam;
/* Computing 2nd power */
	d__3 = q0;
/* Computing 2nd power */
	d__4 = rlam;
	s = log(log(d__1 * d__1 / (d__2 * d__2)) / log(d__3 * d__3 / (d__4 * 
		d__4)));
/* Computing 2nd power */
	d__1 = *qs;
/* Computing 2nd power */
	d__2 = rlam;
/* Computing 2nd power */
	d__3 = rmcc;
/* Computing 2nd power */
	d__4 = rlam;
	sc = log(log(d__1 * d__1 / (d__2 * d__2)) / log(d__3 * d__3 / (d__4 * 
		d__4)));
/* Computing 2nd power */
	d__1 = *qs;
/* Computing 2nd power */
	d__2 = rlam;
/* Computing 2nd power */
	d__3 = rmbb;
/* Computing 2nd power */
	d__4 = rlam;
	sb = log(log(d__1 * d__1 / (d__2 * d__2)) / log(d__3 * d__3 / (d__4 * 
		d__4)));
/* ---------------------- NLO PION ------------------------------ */
	b[0] = 3.73331;
	b[1] = -.74159;
	b[2] = 2.33092;
	a[0] = -3.16946;
	a[1] = -.47683;
	a[2] = .7027;
	a[3] = -.51377;
	a[4] = -.19705;
	a[5] = -.17917;
	a[6] = 2.03394;
	a[7] = -.50764;
	a[8] = -.08565;
	a[9] = .09466;
	a[10] = -.10222;
/* Computing 2nd power */
	d__1 = s;
/* Computing 3rd power */
	d__2 = s;
/* Computing 2nd power */
	d__4 = s;
/* Computing 3rd power */
	d__5 = s;
	d__3 = b[1] + a[3] * s + a[4] * (d__4 * d__4) + a[5] * (d__5 * (d__5 *
		 d__5));
	d__6 = 1. - *x;
/* Computing 2nd power */
	d__8 = s;
/* Computing 3rd power */
	d__9 = s;
	d__7 = b[2] + a[6] * s + a[7] * (d__8 * d__8) + a[8] * (d__9 * (d__9 *
		 d__9));
/* Computing 2nd power */
	d__10 = s;
	dpg = (b[0] + a[0] * s + a[1] * (d__1 * d__1) + a[2] * (d__2 * (d__2 *
		 d__2))) * pow__dd(x, &d__3) * pow__dd(&d__6, &d__7) * ((a[9] * 
		s + a[10] * (d__10 * d__10)) / *x + 1.);
	b[0] = .44809;
	b[1] = -1.47598;
	b[2] = .91338;
	a[0] = -.13828;
	a[1] = -.06951;
	a[2] = .01354;
	a[3] = -.30498;
	a[4] = -.01863;
	a[5] = -.12529;
	a[6] = .64145;
	a[7] = .0727;
	a[8] = -.16989;
	a[9] = .07396;
	a[10] = -.07757;
/* Computing 2nd power */
	d__1 = s;
/* Computing 3rd power */
	d__2 = s;
/* Computing 2nd power */
	d__4 = s;
/* Computing 3rd power */
	d__5 = s;
	d__3 = b[1] + a[3] * s + a[4] * (d__4 * d__4) + a[5] * (d__5 * (d__5 *
		 d__5));
	d__6 = 1. - *x;
/* Computing 2nd power */
	d__8 = s;
/* Computing 3rd power */
	d__9 = s;
	d__7 = b[2] + a[6] * s + a[7] * (d__8 * d__8) + a[8] * (d__9 * (d__9 *
		 d__9));
/* Computing 2nd power */
	d__10 = s;
	dpu = (b[0] + a[0] * s + a[1] * (d__1 * d__1) + a[2] * (d__2 * (d__2 *
		 d__2))) * pow__dd(x, &d__3) * pow__dd(&d__6, &d__7) * ((a[9] * 
		s + a[10] * (d__10 * d__10)) / *x + 1.);
	b[0] = 16.5987;
	b[1] = .13345;
	b[2] = 5.89903;
	a[0] = -18.3856;
	a[1] = 2.44225;
	a[2] = 2.13225;
	a[3] = .22712;
	a[4] = -.83625;
	a[5] = .38526;
	a[6] = -.16911;
	a[7] = .59886;
	a[8] = -.2563;
	a[9] = -.18619;
	a[10] = .87362;
/* Computing 2nd power */
	d__1 = s;
/* Computing 3rd power */
	d__2 = s;
/* Computing 2nd power */
	d__4 = s;
/* Computing 3rd power */
	d__5 = s;
	d__3 = b[1] + a[3] * s + a[4] * (d__4 * d__4) + a[5] * (d__5 * (d__5 *
		 d__5));
	d__6 = 1. - *x;
/* Computing 2nd power */
	d__8 = s;
/* Computing 3rd power */
	d__9 = s;
	d__7 = b[2] + a[6] * s + a[7] * (d__8 * d__8) + a[8] * (d__9 * (d__9 *
		 d__9));
/* Computing 2nd power */
	d__10 = s;
	dps = (b[0] + a[0] * s + a[1] * (d__1 * d__1) + a[2] * (d__2 * (d__2 *
		 d__2))) * pow__dd(x, &d__3) * pow__dd(&d__6, &d__7) * ((a[9] * 
		s + a[10] * (d__10 * d__10)) / *x + 1.);
	b[0] = 6.17173;
	b[1] = -.53618;
	b[2] = 5.60108;
	a[0] = -4.8245;
	a[1] = -1.30844;
	a[2] = 1.95527;
	a[3] = -.27879;
	a[4] = -.51337;
	a[5] = .109;
	a[6] = .83571;
	a[7] = -1.15141;
	a[8] = .77027;
	a[9] = .09268;
	a[10] = -.11267;
/* Computing 2nd power */
	d__1 = sc;
/* Computing 3rd power */
	d__2 = sc;
/* Computing 2nd power */
	d__4 = sc;
/* Computing 3rd power */
	d__5 = sc;
	d__3 = b[1] + a[3] * sc + a[4] * (d__4 * d__4) + a[5] * (d__5 * (d__5 
		* d__5));
	d__6 = 1. - *x;
/* Computing 2nd power */
	d__8 = sc;
/* Computing 3rd power */
	d__9 = sc;
	d__7 = b[2] + a[6] * sc + a[7] * (d__8 * d__8) + a[8] * (d__9 * (d__9 
		* d__9));
/* Computing 2nd power */
	d__10 = sc;
	dpc = (b[0] + a[0] * sc + a[1] * (d__1 * d__1) + a[2] * (d__2 * (d__2 
		* d__2))) * pow__dd(x, &d__3) * pow__dd(&d__6, &d__7) * ((a[9] *
		 sc + a[10] * (d__10 * d__10)) / *x + 1.);
	b[0] = .25944;
	b[1] = -1.98713;
	b[2] = 3.52857;
	a[0] = -.11449;
	a[1] = .03733;
	a[2] = -.18028;
	a[3] = -.35858;
	a[4] = .22277;
	a[5] = -.66413;
	a[6] = .72303;
	a[7] = .4626;
	a[8] = -.99235;
	a[9] = -.02701;
	a[10] = -.02089;
/* Computing 2nd power */
	d__1 = sb;
/* Computing 3rd power */
	d__2 = sb;
/* Computing 2nd power */
	d__4 = sb;
/* Computing 3rd power */
	d__5 = sb;
	d__3 = b[1] + a[3] * sb + a[4] * (d__4 * d__4) + a[5] * (d__5 * (d__5 
		* d__5));
	d__6 = 1. - *x;
/* Computing 2nd power */
	d__8 = sb;
/* Computing 3rd power */
	d__9 = sb;
	d__7 = b[2] + a[6] * sb + a[7] * (d__8 * d__8) + a[8] * (d__9 * (d__9 
		* d__9));
/* Computing 2nd power */
	d__10 = sb;
	dpb = (b[0] + a[0] * sb + a[1] * (d__1 * d__1) + a[2] * (d__2 * (d__2 
		* d__2))) * pow__dd(x, &d__3) * pow__dd(&d__6, &d__7) * ((a[9] *
		 sb + a[10] * (d__10 * d__10)) / *x + 1.);
/* ---------------------- NLO KAON ------------------------------ */
	b[0] = .2314;
	b[1] = -1.364;
	b[2] = 1.79761;
	a[0] = -.33644;
	a[1] = .16204;
	a[2] = -.02598;
	a[3] = .97182;
	a[4] = -.02908;
	a[5] = -.43195;
	a[6] = 1.57116;
	a[7] = .71847;
	a[8] = -.68331;
	a[9] = .36906;
	a[10] = 2.3906;
/* Computing 2nd power */
	d__1 = s;
/* Computing 3rd power */
	d__2 = s;
/* Computing 2nd power */
	d__4 = s;
/* Computing 3rd power */
	d__5 = s;
	d__3 = b[1] + a[3] * s + a[4] * (d__4 * d__4) + a[5] * (d__5 * (d__5 *
		 d__5));
	d__6 = 1. - *x;
/* Computing 2nd power */
	d__8 = s;
/* Computing 3rd power */
	d__9 = s;
	d__7 = b[2] + a[6] * s + a[7] * (d__8 * d__8) + a[8] * (d__9 * (d__9 *
		 d__9));
/* Computing 2nd power */
	d__10 = s;
	dkg = (b[0] + a[0] * s + a[1] * (d__1 * d__1) + a[2] * (d__2 * (d__2 *
		 d__2))) * pow__dd(x, &d__3) * pow__dd(&d__6, &d__7) * ((a[9] * 
		s + a[10] * (d__10 * d__10)) / *x + 1.);
	b[0] = .17806;
	b[1] = -.53733;
	b[2] = .7594;
	a[0] = -.10988;
	a[1] = -.02524;
	a[2] = .03142;
	a[3] = -.60058;
	a[4] = .07863;
	a[5] = .13276;
	a[6] = .61356;
	a[7] = -.43886;
	a[8] = .23942;
	a[9] = .10742;
	a[10] = .128;
/* Computing 2nd power */
	d__1 = s;
/* Computing 3rd power */
	d__2 = s;
/* Computing 2nd power */
	d__4 = s;
/* Computing 3rd power */
	d__5 = s;
	d__3 = b[1] + a[3] * s + a[4] * (d__4 * d__4) + a[5] * (d__5 * (d__5 *
		 d__5));
	d__6 = 1. - *x;
/* Computing 2nd power */
	d__8 = s;
/* Computing 3rd power */
	d__9 = s;
	d__7 = b[2] + a[6] * s + a[7] * (d__8 * d__8) + a[8] * (d__9 * (d__9 *
		 d__9));
/* Computing 2nd power */
	d__10 = s;
	dku = (b[0] + a[0] * s + a[1] * (d__1 * d__1) + a[2] * (d__2 * (d__2 *
		 d__2))) * pow__dd(x, &d__3) * pow__dd(&d__6, &d__7) * ((a[9] * 
		s + a[10] * (d__10 * d__10)) / *x + 1.);
	b[0] = 4.96269;
	b[1] = .05562;
	b[2] = 2.79926;
	a[0] = 1.54098;
	a[1] = -9.06376;
	a[2] = 4.94791;
	a[3] = 1.8866;
	a[4] = -2.9435;
	a[5] = 1.04227;
	a[6] = 3.02991;
	a[7] = -4.14807;
	a[8] = 1.91494;
	a[9] = .8545;
	a[10] = -.61016;
/* Computing 2nd power */
	d__1 = s;
/* Computing 3rd power */
	d__2 = s;
/* Computing 2nd power */
	d__4 = s;
/* Computing 3rd power */
	d__5 = s;
	d__3 = b[1] + a[3] * s + a[4] * (d__4 * d__4) + a[5] * (d__5 * (d__5 *
		 d__5));
	d__6 = 1. - *x;
/* Computing 2nd power */
	d__8 = s;
/* Computing 3rd power */
	d__9 = s;
	d__7 = b[2] + a[6] * s + a[7] * (d__8 * d__8) + a[8] * (d__9 * (d__9 *
		 d__9));
/* Computing 2nd power */
	d__10 = s;
	dkd = (b[0] + a[0] * s + a[1] * (d__1 * d__1) + a[2] * (d__2 * (d__2 *
		 d__2))) * pow__dd(x, &d__3) * pow__dd(&d__6, &d__7) * ((a[9] * 
		s + a[10] * (d__10 * d__10)) / *x + 1.);
	b[0] = 4.25954;
	b[1] = -.24144;
	b[2] = 4.21265;
	a[0] = -5.44309;
	a[1] = 6.11031;
	a[2] = -3.13973;
	a[3] = -1.07757;
	a[4] = 1.52364;
	a[5] = -.74308;
	a[6] = .2559;
	a[7] = .98423;
	a[8] = -.52839;
	a[9] = -.04;
	a[10] = .08695;
/* Computing 2nd power */
	d__1 = sc;
/* Computing 3rd power */
	d__2 = sc;
/* Computing 2nd power */
	d__4 = sc;
/* Computing 3rd power */
	d__5 = sc;
	d__3 = b[1] + a[3] * sc + a[4] * (d__4 * d__4) + a[5] * (d__5 * (d__5 
		* d__5));
	d__6 = 1. - *x;
/* Computing 2nd power */
	d__8 = sc;
/* Computing 3rd power */
	d__9 = sc;
	d__7 = b[2] + a[6] * sc + a[7] * (d__8 * d__8) + a[8] * (d__9 * (d__9 
		* d__9));
/* Computing 2nd power */
	d__10 = sc;
	dkc = (b[0] + a[0] * sc + a[1] * (d__1 * d__1) + a[2] * (d__2 * (d__2 
		* d__2))) * pow__dd(x, &d__3) * pow__dd(&d__6, &d__7) * ((a[9] *
		 sc + a[10] * (d__10 * d__10)) / *x + 1.);
	b[0] = 1.32443;
	b[1] = -.88351;
	b[2] = 6.15221;
	a[0] = -1.41156;
	a[1] = -.04809;
	a[2] = .79066;
	a[3] = -.44818;
	a[4] = -.60073;
	a[5] = .45526;
	a[6] = .46679;
	a[7] = -.50792;
	a[8] = .67006;
	a[9] = -.00477;
	a[10] = -.05503;
/* Computing 2nd power */
	d__1 = sb;
/* Computing 3rd power */
	d__2 = sb;
/* Computing 2nd power */
	d__4 = sb;
/* Computing 3rd power */
	d__5 = sb;
	d__3 = b[1] + a[3] * sb + a[4] * (d__4 * d__4) + a[5] * (d__5 * (d__5 
		* d__5));
	d__6 = 1. - *x;
/* Computing 2nd power */
	d__8 = sb;
/* Computing 3rd power */
	d__9 = sb;
	d__7 = b[2] + a[6] * sb + a[7] * (d__8 * d__8) + a[8] * (d__9 * (d__9 
		* d__9));
/* Computing 2nd power */
	d__10 = sb;
	dkb = (b[0] + a[0] * sb + a[1] * (d__1 * d__1) + a[2] * (d__2 * (d__2 
		* d__2))) * pow__dd(x, &d__3) * pow__dd(&d__6, &d__7) * ((a[9] *
		 sb + a[10] * (d__10 * d__10)) / *x + 1.);
/* ---------------------- NLO PROTON ----------------------------- */
	b[0] = 1.56255;
	b[1] = .01567;
	b[2] = 3.57583;
	a[0] = -1.48158;
	a[1] = -.39439;
	a[2] = .51249;
	a[3] = -2.16232;
	a[4] = 2.47127;
	a[5] = -.93259;
	a[6] = 3.33958;
	a[7] = -3.05265;
	a[8] = 1.21042;
	a[9] = -.84816;
	a[10] = 1.23583;
/* Computing 2nd power */
	d__1 = s;
/* Computing 3rd power */
	d__2 = s;
/* Computing 2nd power */
	d__4 = s;
/* Computing 3rd power */
	d__5 = s;
	d__3 = b[1] + a[3] * s + a[4] * (d__4 * d__4) + a[5] * (d__5 * (d__5 *
		 d__5));
	d__6 = 1. - *x;
/* Computing 2nd power */
	d__8 = s;
/* Computing 3rd power */
	d__9 = s;
	d__7 = b[2] + a[6] * s + a[7] * (d__8 * d__8) + a[8] * (d__9 * (d__9 *
		 d__9));
/* Computing 2nd power */
	d__10 = s;
	dprg = (b[0] + a[0] * s + a[1] * (d__1 * d__1) + a[2] * (d__2 * (d__2 
		* d__2))) * pow__dd(x, &d__3) * pow__dd(&d__6, &d__7) * ((a[9] *
		 s + a[10] * (d__10 * d__10)) / *x + 1.);
	b[0] = 1.25946;
	b[1] = .07124;
	b[2] = 4.12795;
	a[0] = -1.17505;
	a[1] = .3755;
	a[2] = -.01416;
	a[3] = -.29533;
	a[4] = -.2454;
	a[5] = .16543;
	a[6] = .98867;
	a[7] = -.46846;
	a[8] = .2075;
	a[9] = .18957;
	a[10] = -.01116;
/* Computing 2nd power */
	d__1 = s;
/* Computing 3rd power */
	d__2 = s;
/* Computing 2nd power */
	d__4 = s;
/* Computing 3rd power */
	d__5 = s;
	d__3 = b[1] + a[3] * s + a[4] * (d__4 * d__4) + a[5] * (d__5 * (d__5 *
		 d__5));
	d__6 = 1. - *x;
/* Computing 2nd power */
	d__8 = s;
/* Computing 3rd power */
	d__9 = s;
	d__7 = b[2] + a[6] * s + a[7] * (d__8 * d__8) + a[8] * (d__9 * (d__9 *
		 d__9));
/* Computing 2nd power */
	d__10 = s;
	dpru = (b[0] + a[0] * s + a[1] * (d__1 * d__1) + a[2] * (d__2 * (d__2 
		* d__2))) * pow__dd(x, &d__3) * pow__dd(&d__6, &d__7) * ((a[9] *
		 s + a[10] * (d__10 * d__10)) / *x + 1.);
	b[0] = 4.01135;
	b[1] = .17258;
	b[2] = 5.20766;
	a[0] = 8.67124;
	a[1] = -22.7888;
	a[2] = 11.472;
	a[3] = 4.57608;
	a[4] = -9.64835;
	a[5] = 4.61792;
	a[6] = 7.25144;
	a[7] = -12.6313;
	a[8] = 6.07314;
	a[9] = .16931;
	a[10] = -.09541;
/* Computing 2nd power */
	d__1 = s;
/* Computing 3rd power */
	d__2 = s;
/* Computing 2nd power */
	d__4 = s;
/* Computing 3rd power */
	d__5 = s;
	d__3 = b[1] + a[3] * s + a[4] * (d__4 * d__4) + a[5] * (d__5 * (d__5 *
		 d__5));
	d__6 = 1. - *x;
/* Computing 2nd power */
	d__8 = s;
/* Computing 3rd power */
	d__9 = s;
	d__7 = b[2] + a[6] * s + a[7] * (d__8 * d__8) + a[8] * (d__9 * (d__9 *
		 d__9));
/* Computing 2nd power */
	d__10 = s;
	dprs = (b[0] + a[0] * s + a[1] * (d__1 * d__1) + a[2] * (d__2 * (d__2 
		* d__2))) * pow__dd(x, &d__3) * pow__dd(&d__6, &d__7) * ((a[9] *
		 s + a[10] * (d__10 * d__10)) / *x + 1.);
	b[0] = .0825;
	b[1] = -1.6129;
	b[2] = 2.01255;
	a[0] = -.04512;
	a[1] = -.00565;
	a[2] = .009;
	a[3] = -.38012;
	a[4] = -.0684;
	a[5] = .08888;
	a[6] = .63782;
	a[7] = -.14146;
	a[8] = .06083;
	a[9] = -.02958;
	a[10] = .0113;
/* Computing 2nd power */
	d__1 = sc;
/* Computing 3rd power */
	d__2 = sc;
/* Computing 2nd power */
	d__4 = sc;
/* Computing 3rd power */
	d__5 = sc;
	d__3 = b[1] + a[3] * sc + a[4] * (d__4 * d__4) + a[5] * (d__5 * (d__5 
		* d__5));
	d__6 = 1. - *x;
/* Computing 2nd power */
	d__8 = sc;
/* Computing 3rd power */
	d__9 = sc;
	d__7 = b[2] + a[6] * sc + a[7] * (d__8 * d__8) + a[8] * (d__9 * (d__9 
		* d__9));
/* Computing 2nd power */
	d__10 = sc;
	dprc = (b[0] + a[0] * sc + a[1] * (d__1 * d__1) + a[2] * (d__2 * (
		d__2 * d__2))) * pow__dd(x, &d__3) * pow__dd(&d__6, &d__7) * ((
		a[9] * sc + a[10] * (d__10 * d__10)) / *x + 1.);
	b[0] = 24.2916;
	b[1] = .57939;
	b[2] = 12.1207;
	a[0] = -88.3524;
	a[1] = 93.1056;
	a[2] = -17.4089;
	a[3] = -.80783;
	a[4] = -5.072;
	a[5] = -2.45377;
	a[6] = -3.2737;
	a[7] = 1.21188;
	a[8] = -5.50374;
	a[9] = .14628;
	a[10] = -.78634;
/* Computing 2nd power */
	d__1 = sb;
/* Computing 3rd power */
	d__2 = sb;
/* Computing 2nd power */
	d__4 = sb;
/* Computing 3rd power */
	d__5 = sb;
	d__3 = b[1] + a[3] * sb + a[4] * (d__4 * d__4) + a[5] * (d__5 * (d__5 
		* d__5));
	d__6 = 1. - *x;
/* Computing 2nd power */
	d__8 = sb;
/* Computing 3rd power */
	d__9 = sb;
	d__7 = b[2] + a[6] * sb + a[7] * (d__8 * d__8) + a[8] * (d__9 * (d__9 
		* d__9));
/* Computing 2nd power */
	d__10 = sb;
	dprb = (b[0] + a[0] * sb + a[1] * (d__1 * d__1) + a[2] * (d__2 * (
		d__2 * d__2))) * pow__dd(x, &d__3) * pow__dd(&d__6, &d__7) * ((
		a[9] * sb + a[10] * (d__10 * d__10)) / *x + 1.);
    }
/* --- Evaluate different contributions */
    dpd = dpu;
    dks = dku;
    dprd = dpru * .5;
    dk0s = dk0d;
    if (*qs < rmbb) {
	dpb = 0.;
	dkb = 0.;
	dprb = 0.;
	dk0b = 0.;
    }
    if (*qs < rmcc) {
	dpc = 0.;
	dkc = 0.;
	dprc = 0.;
	dk0c = 0.;
    }
    if (*ih == 1) {
	dh[0] = dpg / 2.;
	dh[1] = dpu / 2.;
	dh[2] = dpu / 2.;
	dh[3] = dpd / 2.;
	dh[4] = dpd / 2.;
	dh[5] = dps / 2.;
	dh[6] = dps / 2.;
	dh[7] = dpc / 2.;
	dh[8] = dpc / 2.;
	dh[9] = dpb / 2.;
	dh[10] = dpb / 2.;
    } else if (*ih == 2) {
	dh[0] = dkg / 2.;
	dh[1] = dku / 2.;
	dh[2] = dku / 2.;
	dh[3] = dkd / 2.;
	dh[4] = dkd / 2.;
	dh[5] = dks / 2.;
	dh[6] = dks / 2.;
	dh[7] = dkc / 2.;
	dh[8] = dkc / 2.;
	dh[9] = dkb / 2.;
	dh[10] = dkb / 2.;
    } else if (*ih == 3) {
	dh[0] = dkg / 2.;
	dh[1] = dkd / 2.;
	dh[2] = dkd / 2.;
	dh[3] = dku / 2.;
	dh[4] = dku / 2.;
	dh[5] = dks / 2.;
	dh[6] = dks / 2.;
	dh[7] = dkc / 2.;
	dh[8] = dkc / 2.;
	dh[9] = dkb / 2.;
	dh[10] = dkb / 2.;
    } else if (*ih == 4) {
	dh[0] = dprg / 2.;
	dh[1] = dpru / 2.;
	dh[2] = dpru / 2.;
	dh[3] = dprd / 2.;
	dh[4] = dprd / 2.;
	dh[5] = dprs / 2.;
	dh[6] = dprs / 2.;
	dh[7] = dprc / 2.;
	dh[8] = dprc / 2.;
	dh[9] = dprb / 2.;
	dh[10] = dprb / 2.;
    } else if (*ih == 5) {
	dh[0] = dpg / 2.;
	dh[1] = dpu / 2.;
	dh[2] = dpu / 2.;
	dh[3] = dpd / 2.;
	dh[4] = dpd / 2.;
	dh[5] = dps / 2.;
	dh[6] = dps / 2.;
	dh[7] = dpc / 2.;
	dh[8] = dpc / 2.;
	dh[9] = dpb / 2.;
	dh[10] = dpb / 2.;
    } else if (*ih == 6) {
	dh[0] = dprg / 2.;
	dh[1] = dpru / 4.;
	dh[2] = dpru / 4.;
	dh[3] = dprd;
	dh[4] = dprd;
	dh[5] = dprs / 2.;
	dh[6] = dprs / 2.;
	dh[7] = dprc / 2.;
	dh[8] = dprc / 2.;
	dh[9] = dprb / 2.;
	dh[10] = dprb / 2.;
    } else {
	dh[0] = dpg + dkg + dprg;
	dh[1] = dpu + dku + dpru;
	dh[2] = dpu + dku + dpru;
	dh[3] = dpd + dkd + dprd;
	dh[4] = dpd + dkd + dprd;
	dh[5] = dps + dks + dprs;
	dh[6] = dps + dks + dprs;
	dh[7] = dpc + dkc + dprc;
	dh[8] = dpc + dkc + dprc;
	dh[9] = dpb + dkb + dprb;
	dh[10] = dpb + dkb + dprb;
    }
} /* kkp_ */
