#ifndef Bases_h
#define Bases_h

      //parameter (mxdim = 50, ndmx = 50, leng = 32768)

#include <sstream>

#define mxdim 50
#define ndmx 50
#define leng  32768
#define itm   50

class BSRand;

class Bases
{

private:
    //static const int mxdim,ndmx,leng;
    //static const int itm;

    //common /base0/ jflag,ibases
    int jflag,ibases;

    //common /base1/ xl(mxdim),xu(mxdim),ndim,nwild,ig(mxdim),ncall
    double xl[mxdim],xu[mxdim];
    int ndim,nwild,ig[mxdim],ncall;

    //common /base2/ acc1,acc2,itmx1,itmx2
    double acc1,acc2;
    int   itmx1,itmx2;

    //common /base3/ scalls,wgt,ti,tsi,tacc,it
    double scalls,wgt,ti,tsi,tacc;
    int it;

    //common /base4/ xi(ndmx,mxdim),dx(mxdim),dxd(leng),dxp(leng),
    // .               nd,ng,npg,ma(mxdim)
    double xi[ndmx][mxdim],dx[mxdim],dxd[leng],dxp[leng];
    int nd,ng,npg,ma[mxdim];
    // xi[nd][ndim]

    //real*4 time, eff, wrong, trslt, tstd, pcnt
    //common /base5/ itrat(itm,0:1),time(itm,0:2),eff(itm,0:1),
    //               wrong(itm,0:1),reslt(itm,0:1),acstd(itm,0:1),
    //               trslt(itm,0:1),tstd(itm,0:1),pcnt(itm,0:1)
    int    itrat[itm][2];
    float  time[itm][3],eff[itm][2],wrong[itm][2];
    double reslt[itm][2],acstd[itm][2];
    float  trslt[itm][2],tstd[itm][2],pcnt[itm][2];

    //common /base6/ d(ndmx,mxdim),alph,xsave(ndmx,mxdim),xti,xtsi,xacc,itsx
    double d[ndmx][mxdim],alph,xsave[ndmx][mxdim],xti,xtsi,xacc;
    int itsx; // <=igopt

    //real*4 stime
    //common /bsrslt/avgi,sd,chi2a,stime,itg,itf
    //common /bsrslt/avgi,sd,chi2a,stime,it1,itf
    double avgi,sd,chi2a;
    float  stime;
    int    itg,itf;

    //     intv = ( 0 / 1 / any ) = ( Batch / Batch(Unix) / Interactive )
    //     ipnt = ( 0 / any ) = ( IBM Type / Ascii printer )
    //common /bscntl/ intv, ipnt, nloop, mloop
    int intv, ipnt, nloop, mloop;


      //common /xhcntl/ lock
    int lock;

    //common/ninfo/ nodeid, numnod
    int nodeid, numnod;

    //real*4 timebs,timint,timesp,time0,rtime,timeb1,timeb2,times1
    //common /btime1/ time0,rtime,timeb1,timeb2,times1
    float time0,rtime,timeb1,timeb2,times1;


    //common /btime2/ timebs(0:2),timint,timesp(0:2)
    float timebs[3],timint,timesp[3];

     //common /bdate/ idate(3),itime(2)
    //            IDATE(1) : year        ITIME(1) : hour
    //            IDATE(2) : month       ITIME(2) : minute
    //            IDATE(3) : day
    int idate[3],itime[2];


    //character*80 error
    //common /bwarn1/ nerror
    //common /bwarn2/ error(3,3)
    int nerror;
    std::ostringstream error[3][3];

    BSRand* bsrand;

public:
    Bases();
    virtual ~Bases();
    double bases();
    void bsinit();
    virtual double func(double* x)=0;
    double getError() {return sd;}
    void setPrint(int i) {intv = i; }
    void setXL(int i, double a) { xl[i]=a;}
    void setXU(int i, double a) { xu[i]=a;}
    void setIG(int i, int a) { ig[i]=a;}
    void setNdim(int i) {ndim=i;}
    void setNoOfSample(int i) {ncall=i;}
    void setNWild(int i) {nwild=i;}
    void setIteration1(double a, int i) {acc1=a; itmx1=i;}
    void setIteration2(double a, int i) {acc2=a; itmx2=i;}
    void setTuneValue(double a) {alph = a;}
    void defineVariable(int n, int nw, double* a, double* b,int* c) {
	ndim = n;
	nwild=nw;
	for(int i=0;i<ndim;i++) {
	    xl[i]=a[i];
	    xu[i]=b[i];
	    ig[i]=c[i];
	}
    }


private:
    void bschck();
    void bsdate();
    void bsetgu();
    void bsintg();
    void bsprnt( int id, int ip1=0, int ip2=0 );
    void bslist(int i, int istep );
    void bsutim( int job, int id );
    void bhrset();
    void bsetgv( int i );
    void bsordr( double are, double& fx2, double& order, int& iordr );
    void bstcnv( float time, int& ih, int& mn, int& is1, int& is2 );
    //void bhsave();

};

//#ifndef BSRand_h
//#define BSRand_h

class BSRand
{
private:
    int seed;
    float rdm[31],rm1,rm2;
    int ia1,ic1,m1,ix1;
    int ia2,ic2,m2,ix2;
    int ia3,ic3,m3,ix3;

public:
    void drnset( int iseed );
    double drn();
};

//#endif

class BSTime
{
private:
    static float time_init;

public:
    static float bstime(int iflg );
    static void bsdate(int* idate, int* itime);

private:
    static float uxtime();
    static void uxdate(int& year,int& mon,int& day,int& hour,int& min);
};

#endif
