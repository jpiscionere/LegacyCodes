#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <err.h>
#include "/home/piscioja/lib/rng_gsl.c"


#define sqr(x) ((x)*(x))
#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))
#define mabs(A) ((A) < 0.0 ? -(A) : (A))
#define cnint(x) ((x-floor(x)) < 0.5 ? floor(x) : ceil(x))
#define csign(x) (x < 0.0 ? -1 : 1)
#define omega_m (0.25)
#define G (4.32E-9) // (km/s)^2 Mpc/Msun
#define Ho (100.0)
#define omega_l (0.75)
#define PI (4.0 * atan(1.0))
/* this is the maximum number of bins in the NFW profile calculations */
#define NMAX_NFW 10000
/* multiplicative factor on predicted number of galaxies wrt Nhalos */
#define GALAXY_COUNT_FACTOR 1.5

/* just declaring this as a file global, since it's read in.  Shouldn't really do this.. */
static double M_STAR;

void Printhelp(void)
{
  fprintf(stderr, "%s", "\n"
  "  --- halobias_so_nfw Ncenopt Nsatopt PNNopt logMmin siglogM logM0 logM1 alpha center gamma fgal simulation Deltavir M_STAR redshift PNMfile seed <.halos > Galaxies\n"
  "  --- Creates a biased galaxy distribution by populating a N-body halo catalog using an input HOD.\n"
  "  --- This version does not read particle info, but places galaxies around halo centers according to an adopted profile.\n"
  "  --- The input halo catolgue is in ascii format.\n"
  "\n"
  "     * Ncenopt = 0 : Ncen = 0\n"
  "               = 1 : Ncen = 1 for M>Mmin\n"
  "               = 2 : <Ncen> = exp(-Mmin/M)                                               (Zehavi et al. 2005; Tinker et al. 2005)\n"
  "               = 3 : <Ncen> = 0.5*[1 + erf((logM-logMmin)/siglogM)]                      (Zheng, Coil, & Zehavi 2007)\n"
  "\n"
  "     * Nsatopt = 0 : Nsat = 0 ;\n"
  "               = 1 : <Nsat> = (M/M1)^alpha for M>Mmin                                    (Kravtsov et al. 2004)\n"
  "               = 2 : <Nsat> = (M/M1)^alpha for M>M0\n"
  "               = 3 : <Nsat> = exp(-M0/(M-Mmin))*(M/M1)^alpha                             (Zehavi et al. 2005; Tinker et al. 2005)\n"
  "               = 4 : <Nsat> = 0.5*[1 + erf((logM-logMmin)/siglogM)] * ((M-M0)/M1)^alpha  (Zheng, Coil & Zehavi 2007)\n"
  "\n"
  "     * Pnnopt  = 0 : Average   (Nsat = nint(<Nsat>) - with frequencies to preserve mean)\n"
  "               = 1 : Poisson   (Nsat drawn from Poisson distribution)\n"
  "               = 2 : Binomial  (Nsat drawn from a Binomial distribution)\n"
  "               = 3 : Negative Binomial (Nsat drawn from a Negative Binomial distribution)\n"
  "\n"
  "     * logMmin = minimum mass of halo that can contain a galaxy (in units of Msun/h)\n"
  "     * siglogM = width of cutoff at Mmin (scatter in the mass-luminosity relation)\n"
  "     * logM0   = minimum halo mass that can contain a satellite galaxy (in units of Msun/h)\n"
  "     * logM1   = halo mass for which <Nsat>=1 (in units of Msun/h)\n"
  "     * alpha   = slope of the <Nsat> - M relation\n"
  "\n"
  "     * center  = 0 : Don't force the \"central\" galaxy to be at center of its halo\n"
  "               = 1 : Force the \"central\" galaxy to be at center of its halo\n"
  "\n"
  "     * gamma   = What gamma to use in the profile (default = 1)\n"
  "     * fgal    = Difference between dark matter and galaxy concentration (default = 1)\n"
  "\n"
  "     * simulation  = The simulation the .halos file is from.  Gives the values for mass of particles, Mstar, size of box & redshift of box.\n"
  "                   = -18.0  : Consuelo < -18\n"
  "                   = -19.0  : Consuelo < -19\n"
  "                   = -20.0  : Esmeralda < -20\n"
  "                   = -21.0  : Carmen < -21\n"
  "\n"
  "     * Deltavir= Virial overdensity.  Must be consistent with fof linking length!\n"
  "     * M_Star  = Value of Mstar mass. \n"
  "               Median Redshift values:\n"
  "               = 2.49e12 : Consuelo\n"
  "               = 2.29e12 : Esmeralda\n"
  "               = 1.97e12 : Carmen\n"
  "     * redshift=  Redshift of snapshots.  Should be consistent with M_star. \n"
  "                 Median Redshift values:\n"
  "               = 0.054 : Consuelo\n"
  "               = 0.082 : Esmeralda\n"
  "               = 1.32  : Carmen\n"  
  "\n"
  "     * seed    = seed for random number generator (unsigned long)\n"
  "     * < .halos = the *.halos ascii halo center file containing halo information\n"
  "     * > Biased galaxy distribution (fast food format)\n"
  );
}

/* A workspace for constants, and allocated arrays
 * to aid efficiency for the NFW_radius routine */
typedef struct {
  int nmax;

  void * rng_rad;
  void * rng_pos;

  double * spd; /* allocated array for sum_probability_distribution */
  double * radius;
} NFW_WORK;

NFW_WORK init_nfw_workspace( const unsigned long seed_pos, const unsigned long seed_rad, const int nmax )
{
  NFW_WORK work;

  work.nmax = nmax;

  /* create two independent streams of random numbers,
   * for simplicity: initialize them with same seed
   * (they are still independent) */
  work.rng_rad = rng_create(seed_rad);
  work.rng_pos = rng_create(seed_pos);

  work.spd = (double *) calloc(nmax, sizeof(double));
  assert(work.spd != NULL);
  work.radius = (double *) calloc(nmax, sizeof(double));
  assert(work.radius != NULL);

  return work;
}

void free_nfw_workspace( NFW_WORK * work )
{
  work->nmax = 0;

  rng_free(work->rng_rad);
  rng_free(work->rng_pos);

  assert(work->spd != NULL) ;
  free(work->spd);
  work->spd = NULL;

  assert(work->radius != NULL) ;
  free(work->radius);
  work->radius= NULL;
}

double nfw_density( const double conc, const double R_vir, const double gamma, const double r )
{

  double rho_nfw;
  rho_nfw = 1./(pow(1.0/conc + r / R_vir,(3.0-gamma)) * pow(r/R_vir,gamma)) ;
  return rho_nfw;
}

void NFW_radius(NFW_WORK work, const double M,const double a, const double redshift,
                const double gamma, const double fgal,
                const int nsat, const double R_vir )
{

  int    i,i_nfw,j_nfw,k_nfw;
  double delta_r;
  double conc, rand_nfw;
  double r_0,rho_nfw,probability_distribution,*sum_probability_distribution;
  double radius_nfw;
  double * sat_position;
  const double prefac=4.0/3.0*PI;
  sum_probability_distribution = work.spd;
  sat_position = work.radius;

  // Reset all values, just in case
  for(i=0; i<nsat; i++) {
      sat_position[i] = 0.0;
  }

  //Concentration of halo
  conc = fgal/(1.0+redshift) * 11.0*pow((M/M_STAR),-0.13);

  delta_r = 1.0E-3; // steps of 1 kpc

  r_0 = 1.0e-3; // starting at 1 kpc

  i_nfw=floor((a*R_vir-r_0)/delta_r); // number of bins
  assert(i_nfw < work.nmax);

  // Probability at inner point
  rho_nfw = nfw_density( conc, R_vir, gamma, r_0);
  probability_distribution = 4.0*PI * sqr(r_0) * rho_nfw * delta_r;

  sum_probability_distribution[0] = probability_distribution;
  sum_probability_distribution[i_nfw] = 0 ;

  // Calulate probabilities in radial bins
  for(j_nfw=1;j_nfw<i_nfw;j_nfw++)
    {
      double r_curr = r_0 + delta_r * j_nfw;
      double r_prev = r_0 + delta_r * (j_nfw - 1);
      rho_nfw = nfw_density( conc, R_vir, gamma, r_curr);
      probability_distribution = prefac*( r_curr*r_curr*r_curr - r_prev*r_prev*r_prev) * rho_nfw; 

      sum_probability_distribution[j_nfw] =  sum_probability_distribution[j_nfw - 1] + probability_distribution;
    }

  // Normalize cumulative probability distribution
  for(k_nfw=0;k_nfw<i_nfw;k_nfw++)
    {
      sum_probability_distribution[k_nfw] /= sum_probability_distribution[i_nfw - 1];
    }

  for(i=0;i<nsat;i++)
    {
      rand_nfw = rng_uniform( work.rng_rad );

      k_nfw=0;
      while(rand_nfw > sum_probability_distribution[k_nfw])
        {
          k_nfw++;
          /* this means we're going beyond what we populated, so don't */
          assert(k_nfw < i_nfw);
        }

      radius_nfw = r_0 + delta_r * k_nfw;
      sat_position[i] = radius_nfw;
    }

}

int main(int argc, char *argv[])
{
  int i,j,k ;
/*---Arguments-------------------------*/
  int iarg,Ncenopt,Nsatopt,PNNopt,center ;
  unsigned long seed;
  double logMmin,siglogM,logM0,logM1,alpha,Scale ;
  double Deltavir,gamma,fgal ; /* these are input parameters */
  double Mmin,M0,M1,simulation ;
  FILE *fp2 ;
/*---Halo-centers----------------------*/
  int Nhalos;
  double mp,Lbox,abox,redshift;
  double *Rvir,*Mhalo,*xcen,*ycen,*zcen,*vxcen,*vycen,*vzcen ;
/*---Halo-distribution-variables-------*/
  int halo;
  double Mh,logMh;
/*---Bias-parameters-------------------*/
  int Ncen,Nsat;
  double radius,Theta,n1,Ncenavg,Nsatavg;
/*---Galaxy-distribution-variables-----*/
  int gal,Ngalmax,Ngal,*idat ;
  double *xg,*yg,*zg,*vxg,*vyg,*vzg ;
  float *fdat,znow ;
/*---Functions-------------------------*/
  int Average(double,void *) ;
  int Poisson(double,void *) ;
  int Binomial(double,void *) ;
  int NegBinomial(double,void *) ;
  char * pnmfile;
  NFW_WORK work;
  void * rng;

/*---Read-Arguments-------------------------------------------------------------*/

  if(argc<17)
    {
      Printhelp() ;
      return -1 ;
    }
  sscanf(argv[1],"%d",&Ncenopt) ;
  sscanf(argv[2],"%d",&Nsatopt) ;
  sscanf(argv[3],"%d",&PNNopt) ;
  sscanf(argv[4],"%lf",&logMmin) ;
  Mmin=pow(10.,logMmin) ;
  sscanf(argv[5],"%lf",&siglogM) ;
  sscanf(argv[6],"%lf",&logM0) ;
  M0=pow(10.,logM0) ;
  sscanf(argv[7],"%lf",&logM1) ;
  M1=pow(10.,logM1) ;
  sscanf(argv[8],"%lf",&alpha) ;
  sscanf(argv[9],"%d",&center) ;
  sscanf(argv[10],"%lf",&gamma) ;
  sscanf(argv[11],"%lf",&fgal) ;
  sscanf(argv[12],"%lf",&simulation) ;

{
  if((simulation == -18.0 ) | (simulation == -19.0 ))
  {
    mp=0.187 *1E10;
    Lbox=420.;
    Nhalos=6E6;
  }else if(simulation == -20.0){
    mp=0.931 *1e10;
    Lbox=640.;
    Nhalos=5E6;
  }else if(simulation == -21.0){
    mp=4.938 * 1E10;
    Lbox=1000.;
    Nhalos=4E6;
  }else{
  fprintf(stderr,"halobias> Bad value for simulation! Only -18,-19,-20 & -21 supported!\n");
  return(-1)  ;
  } 


}
  sscanf(argv[13],"%lf",&Deltavir); 
  sscanf(argv[14],"%lf",&M_STAR);
  sscanf(argv[15],"%lf",&redshift);
    abox=1.0/(1.0+redshift); 
 
  sscanf(argv[16],"%lu",&seed) ; // This seed is for choosing number of galaxies

  {
    unsigned long seed1 = 1234; // This seed is for choosing radial positions
    unsigned long seed2 = 9876543210; // This seed is for choosing angular positions
    work = init_nfw_workspace( seed1, seed2, NMAX_NFW );
  }
  rng = rng_create(seed);

  iarg=17 ;

 

/*---Test Arguments ---*/
  if(logMmin <= 0 || logM0 <= 0 || logM1 <= 0) {
    fprintf(stderr,"halobias> Bad values for logMass!  (logMmin, logM0, logM1) = (%g,%g,%g)\n", logMmin, logM0, logM1) ;
    return(-1) ;
  }

// Scale  = Rmax(halo)/Rvir(halo)
  Scale = 1. ;
//  fprintf(stderr,"halobias> The redshift of this snapshot = %lf\n",redshift);
  

//Need to see how long this will take, might have to think of a better way of doing this
  int junki, Npart;
  double junkf;
  double junklf,rho_mean;
  i = j = 0 ;
  rho_mean = 3.0/(8.0*PI)*sqr(Ho)/G * omega_m * pow((1.0 + redshift),3); 
  
/*---Read-group-info-in-BGC-files-----------------------------------------------------*/

  /* we're allocating memory for ALL halos in .halos files */
  Mhalo=(double *) calloc(Nhalos,sizeof(double)) ;
  Rvir=(double *) calloc(Nhalos,sizeof(double)) ;
  xcen=(double *) calloc(Nhalos,sizeof(double)) ;
  ycen=(double *) calloc(Nhalos,sizeof(double)) ;
  zcen=(double *) calloc(Nhalos,sizeof(double)) ;
  vxcen=(double *) calloc(Nhalos,sizeof(double)) ;
  vycen=(double *) calloc(Nhalos,sizeof(double)) ;
  vzcen=(double *) calloc(Nhalos,sizeof(double)) ;
 
  double mass_counter=0; 
  
  i=0; 
  while(fscanf(stdin,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %d",&junki,&Npart,&junklf,&xcen[i],&ycen[i],&zcen[i],&vxcen[i],&vycen[i],&vzcen[i],&junklf,&junki)!=EOF) 
  {

/*---Apply-Warren-et-al-correction-to-mass---*/
    
      Mhalo[i] = mp*((double)Npart*(1.-pow(((double)Npart),-0.6))) ;
      mass_counter += Mhalo[i];
      Rvir[i]= pow(3.0*Mhalo[i]/(4.0*PI*Deltavir*rho_mean),0.33333); 
      if(Mhalo[i]>=Mmin)   
          j++ ;
       
      i++;
    }

fprintf(stderr,"Mass Counter=%e,Mass of Last Halos = %lf, Rvir of last halo= %lf\n",mass_counter,Mhalo[i-1],Rvir[i-1]);
Nhalos=i;



/*---Setup-galaxy-arrays--------------------------------------------------------*/

  Ngalmax = Nhalos * GALAXY_COUNT_FACTOR;
  xg =(double *) calloc(Ngalmax,sizeof(double)) ;
  yg =(double *) calloc(Ngalmax,sizeof(double)) ;
  zg =(double *) calloc(Ngalmax,sizeof(double)) ;
  vxg=(double *) calloc(Ngalmax,sizeof(double)) ;
  vyg=(double *) calloc(Ngalmax,sizeof(double)) ;
  vzg=(double *) calloc(Ngalmax,sizeof(double)) ;

  
 







  gal=0 ;

/*---Loop-over-halos------------------------------------------------------------*/
  for(i=0;i<Nhalos;i++)
    {

      assert(gal < Ngalmax); /* you're trying to make more galaxies than we allocated space for!
				increase GALAXY_COUNT_FACTOR to fix this */
      halo = i;
      Mh = Mhalo[halo] ;
      logMh = log10(Mh) ;

/*---Navg(M)--------------------------------------------------------------------*/

      if(Ncenopt==0)
	{
	  Ncen = 0 ;
	}
      else if(Ncenopt==1)
	{
	  if(Mh<Mmin)
	    Ncen = 0 ;
	  else
	    Ncen = 1 ;
	}
      else if(Ncenopt==2)
	{
	  Ncenavg = exp(-Mmin/Mh) ;
	  Ncen = Average(Ncenavg, rng ) ;
	}
      else if(Ncenopt==3)
	{
	  Ncenavg = 0.5*(1 + erff((logMh-logMmin)/siglogM)) ;
	  Ncen = Average(Ncenavg, rng) ;
	}
      else
	{
	  fprintf(stderr,"halobias> value for Ncenopt is out of range.\n") ;
	  return(-1) ;
	}

/*----------------*/

      if(Ncen==0)
	{
	  Nsatavg = 0. ;
	  Nsat = 0 ;
	}
      else
	{
	  if(Nsatopt==0)
	    {
	      Nsatavg = 0. ;
	      Nsat = 0 ;
	    }
	  else if(Nsatopt==1)
	    {
	      if(Mh<Mmin)
		{
		  Nsatavg = 0. ;
		  Nsat = 0 ;
		}
	      else
		{
		  Nsatavg = pow((Mh/M1),alpha) ;
		  if(PNNopt==0)
		    Nsat = Average(Nsatavg, rng) ;
		  else if(PNNopt==1)
		    Nsat = Poisson(Nsatavg, rng) ;
		  else if(PNNopt==2)
		    Nsat = Binomial(Nsatavg, rng) ;
		  else if(PNNopt==3)
		    Nsat = NegBinomial(Nsatavg, rng) ;
		  else
		    {
		      fprintf(stderr,"halobias> value for PNNopt is out of range.\n") ;
		      return(-1) ;
		    }
		}
	    }
	  else if(Nsatopt==2)
	    {
	      if(Mh<M0)
		{
		  Nsatavg = 0. ;
		  Nsat = 0 ;
		}
	      else
		{
		  Nsatavg = pow((Mh/M1),alpha) ;
		  if(PNNopt==0)
		    Nsat = Average(Nsatavg, rng) ;
		  else if(PNNopt==1)
		    Nsat = Poisson(Nsatavg, rng) ;
		  else if(PNNopt==2)
		    Nsat = Binomial(Nsatavg, rng) ;
		  else if(PNNopt==3)
		    Nsat = NegBinomial(Nsatavg, rng) ;
		  else
		    {
		      fprintf(stderr,"halobias> value for PNNopt is out of range.\n") ;
		      return(-1) ;
		    }
                }
	    }
	  else if(Nsatopt==3)
	    {
	      Nsatavg = exp(-M0/(Mh-Mmin))*pow((Mh/M1),alpha) ;
	      if(PNNopt==0)
		Nsat = Average(Nsatavg, rng) ;
	      else if(PNNopt==1)
		Nsat = Poisson(Nsatavg, rng) ;
	      else if(PNNopt==2)
		Nsat = Binomial(Nsatavg, rng) ;
	      else if(PNNopt==3)
		Nsat = NegBinomial(Nsatavg, rng) ;
	      else
		{
		  fprintf(stderr,"halobias> value for PNNopt is out of range.\n") ;
		  return(-1) ;
		}
	    }
	  else if(Nsatopt==4)
	    {
	      if( Mh > M0 )
		{
		  Nsatavg = 0.5*(1 + erff((logMh-logMmin)/siglogM))*pow(((Mh-M0)/M1),alpha) ;
		  if(PNNopt==0)
		    Nsat = Average(Nsatavg, rng) ;
		  else if(PNNopt==1)
		    Nsat = Poisson(Nsatavg, rng) ;
		  else if(PNNopt==2)
		    Nsat = Binomial(Nsatavg, rng) ;
		  else if(PNNopt==3)
		    Nsat = NegBinomial(Nsatavg, rng) ;
		  else
		    {
		      fprintf(stderr,"halobias> value for PNNopt is out of range.\n") ;
		      return(-1) ;
		    }
		}

	      else
		{
		  Nsat = 0 ;
		}
	    }
	  else
	    {
	      fprintf(stderr,"halobias> value for Nsatopt is out of range.\n") ;
	      return(-1) ;
	    }

	}

/*---Print-out-M-vs-N-in-PNMfile------------------------------------------------*/





/*---Force-a-galaxy-at-halo-center----------------------------------------------*/

      if(Ncen>0)
	{
	  if(center==1)
	    {
	      xg[gal] = xcen[halo] ;
	      yg[gal] = ycen[halo] ;
	      zg[gal] = zcen[halo] ;
	      vxg[gal] = vxcen[halo]*sqrt(abox) ;
	      vyg[gal] = (vycen[halo]*sqrt(abox)) ;
	      vzg[gal] = (vzcen[halo]*sqrt(abox)) ;
	      gal++ ;
	    }
/*---otherwise-deal-with-cen-as-if-it-were-a-sat---*/
	  else
	    {
	      Nsat++ ;
	    }

            //Closes the Ncen >0 loop
	}

/*---Only-continue-if-at-least-one-satellite------------------------------------*/

      if(Nsat>0)
	{
    
/*---Assign Satallites a position---------------------------*/
	  double *sat_position;
	  /* the NFW_radius puts Nsat entries into the work.radius array based on NFW profile */
	  NFW_radius(work, Mh, Scale, redshift, gamma, fgal, Nsat, Rvir[halo]) ;
	  sat_position = work.radius ;

	  for(k=0;k<Nsat;k++)
	    {
	      radius = sat_position[k];
	      n1 = 2 * ( rng_uniform(work.rng_pos) - 0.5 );
	      Theta = 2 * PI * rng_uniform(work.rng_pos);

	      xg[gal] = (double)(radius * ( cos(Theta)*sqrt(1-sqr(n1)) ) + xcen[halo]);

	      if( xg[gal] > (Lbox) )
		{
		  xg[gal] = xg[gal] - Lbox;
		}
	      else if( xg[gal] < 0 )
		{
		  xg[gal] = xg[gal] + Lbox;
		}

	      yg[gal] = (double)(radius * ( sin(Theta)*sqrt(1-sqr(n1)) ) + ycen[halo]);

	      if( yg[gal] > Lbox )
		{
		  yg[gal] = yg[gal] - Lbox;
		}
	      else if( yg[gal] < 0 )
		{
		  yg[gal] = yg[gal] + Lbox;
		}

	      zg[gal] = (double)(n1*radius + zcen[halo]);

	      if( zg[gal] > Lbox )
		{
		  zg[gal] = zg[gal] - Lbox;
		}
	      else if( zg[gal] < 0 )
		{
		  zg[gal] = zg[gal] + Lbox;
		}
        vxg[gal] = (vxcen[halo]*sqrt(abox)) ;
        vyg[gal] = (vycen[halo]*sqrt(abox)) ;
        vzg[gal] = (vzcen[halo]*sqrt(abox)) ;
	      gal++;
	    }

	}
    }
  free_nfw_workspace(&work);
//fprintf(stderr,"WRITE GALAXY FILE \n");
/*---Write-galaxy-file----------------------------------------------------------*/
  

  Ngal = gal ;
  fprintf(stderr,"Ngal = %d\n",Ngal);
//  fprintf(stderr,"halobias> Writing galaxy file. Ngal = %d\n",Ngal) ;

  idat=(int *)calloc(5,sizeof(int)) ;
  fdat=(float *)calloc(9,sizeof(float)) ;

  idat[0] = (int)Lbox ;
  idat[1] = Ngal ;
  idat[2]=idat[3]=idat[4]=0 ;
  fdat[0] = Lbox ;
  fdat[1]=fdat[2]=fdat[3]=fdat[4]=fdat[5]=fdat[6]=fdat[7]=fdat[8]=0.0 ;
  znow=0.0 ;

  ftwrite(idat,sizeof(int),5,stdout);
  ftwrite(fdat,sizeof(float),9,stdout);
  ftwrite(&znow,sizeof(float),1,stdout);
  ftwrite(xg,sizeof(*xg),Ngal,stdout);
  ftwrite(yg,sizeof(*xg),Ngal,stdout);
  ftwrite(zg,sizeof(*xg),Ngal,stdout);
  ftwrite(vxg,sizeof(*xg),Ngal,stdout);
  ftwrite(vyg,sizeof(*xg),Ngal,stdout);
  ftwrite(vzg,sizeof(*xg),Ngal,stdout);

  return 0 ;
}
/******************/
/*   FUNCTIONS    */
/******************/

/*-Average-distribution------------------------------*/
int Average(double Nexp, void * rng )
{
  int Nact ;
  double rand ;

  rand = rng_uniform(rng);

  if(rand<=(Nexp-(int)(Nexp))) Nact = (int)(Nexp+1) ;
  else Nact = (int)(Nexp) ;

  return Nact ;
}

/*-Poisson-distribution------------------------------*/
int Poisson(double Nexp, void * rng )
{
  int i,Nmax,Nact=0 ;
  double x,sigma,P,*Sum,rand ;

  sigma = sqrt(Nexp) ;
  if(Nexp>=0.6) Nmax = (int)(10*sigma+Nexp) ;
  else Nmax = 8;
  Sum=(double *)calloc(Nmax,sizeof(double)) ;

  P = exp(-Nexp) ;
  Sum[0] = P ;
  for(i=1;i<Nmax;i++)
    {
      x = (double)i ;
      P *= Nexp/x ;
      Sum[i] = Sum[i-1]+P ;
    }

  rand = rng_uniform(rng);
  for(i=0;i<Nmax;i++)
    if(rand<=Sum[i])
      {
        Nact = i ;
        break ;
      }
  free(Sum) ;
  return Nact ;
}

/*-Binomial-distribution-----------------------------*/
int Binomial(double Nexp, void * rng)
{
  int i,Nmax,Nact=0,r ;
  double x,p,q,Prob,var,*Sum ;
  double rand ;

  r=(int)(2*Nexp+1) ;
  p=Nexp/(double)r ;
  q=1-p ;
  var=sqrt(r*p*q) ;
  if(Nexp>=4) Nmax = r;
  else Nmax=8 ;
  Sum=(double *)calloc(Nmax,sizeof(double)) ;
  Prob = pow(q,(double)r) ;
  Sum[0] = Prob ;

  for(i=1;i<Nmax;i++)
    {
      x = (double)i ;
      Prob *= ((r-x+1)/x)*p/q ;
      Sum[i] = Sum[i-1]+Prob ;
    }

  rand = rng_uniform(rng);

  for(i=0;i<Nmax;i++)
    if(rand<=Sum[i])
      {
        Nact = i ;
        break ;
      }
  free(Sum) ;
  return Nact ;
}

/*-Negative-Binomial-distribution--------------------*/
int NegBinomial(double Nexp, void * rng )
{
  int i,Nmax,Nact=0,r ;
  double x,p,q,Prob,var,*Sum ;
  double rand ;

  r=(int)Nexp+1 ;
  p=1/(1+Nexp/(double)r) ;
  q=1-p ;
  var=sqrt(r*q/(p*p)) ;

  if(Nexp>=3) Nmax = (int)(10*var+Nexp) ;
  else Nmax = 30;
  Sum=(double *)calloc(Nmax,sizeof(double)) ;

  Prob = pow(p,(double)r) ;
  Sum[0] = Prob ;
  for(i=1;i<Nmax;i++)
    {
      x = (double)i ;
      Prob *= ((r+x-1)/x)*pow(q,x)/pow(q,x-1) ;
      if(i>r)
        {
          if(Prob>1e-20) Prob=Prob ;
          else Prob=0 ;
        }
      Sum[i] = Sum[i-1]+Prob ;
    }

  rand = rng_uniform( rng );
  for(i=0;i<Nmax;i++)
    if(rand<=Sum[i])
      {
        Nact = i ;
        break ;
      }
  free(Sum) ;
  return Nact ;
}

// vim: ts=2 sts=2 sw=2 expandtab
