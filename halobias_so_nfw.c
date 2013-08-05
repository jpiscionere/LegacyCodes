#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <err.h>
#include </home/berlinaa/Source/Lib/nr.h>
#include "/home/berlinaa/Source/BGC_utils/bgc2_read_utils.c"
#include "/home/berlinaa/Source/Lib/rng_gsl.c"

#define sqr(x) ((x)*(x))
#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))
#define mabs(A) ((A) < 0.0 ? -(A) : (A))
#define cnint(x) ((x-floor(x)) < 0.5 ? floor(x) : ceil(x))
#define csign(x) (x < 0.0 ? -1 : 1)

#define PI (4.0 * atan(1.0))
/* this is the maximum number of bins in the NFW profile calculations */
#define NMAX_NFW 10000
/* multiplicative factor on predicted number of galaxies wrt Nhalos */
#define GALAXY_COUNT_FACTOR 1.0

/* just declaring this as a file global, since it's read in.  Shouldn't really do this.. */
static double M_STAR;

void Printhelp(void)
{
  fprintf(stderr, "%s", "\n"
  "  --- halobias_so_nfw Ncenopt Nsatopt PNNopt logMmin siglogM logM0 logM1 alpha center gamma fgal Mstar PNMfile seed *BGCfiles > Galaxies\n"
  "  --- Creates a biased galaxy distribution by populating a N-body halo catalog using an input HOD.\n"
  "  --- This version does not read particle info, but places galaxies around halo centers according to an adopted profile.\n"
  "  --- The input halo catalog consists of SO halos in the BGC2 format.\n"
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
  "     * Mstar   = Value of Mstar mass (non-log10)\n"
  "\n"
  "     * PNMfile = output file containing number of galaxies selected for each halo\n"
  "     * seed    = seed for random number generator (unsigned long)\n"
  "     * BGCfiles = *.bgc2 halo files (with or without PARTICLE_DATA)\n"
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
      probability_distribution = 4.0/3.0*PI*( pow(r_curr,3.0) - pow(r_prev,3.0 )) * rho_nfw;

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
  double Mmin,M0,M1 ;
  FILE *fp1,*fp2 ;
/*---Halo-centers----------------------*/
  int Nhalos;
  double mp,Lbox,abox;
  float *Rvir,*Mhalo,*xcen,*ycen,*zcen,*vxcen,*vycen,*vzcen,redshift ;
/*---Halo-distribution-variables-------*/
  int halo;
  char *bgc_file ;
  OUTPUT_HEADER hdr ;
  double Mh,logMh;
/*---Bias-parameters-------------------*/
  int Ncen,Nsat;
  double radius,Theta,n1,Ncenavg,Nsatavg;
/*---Galaxy-distribution-variables-----*/
  int gal,Ngalmax,Ngal,*idat ;
  float *xg,*yg,*zg,*vxg,*vyg,*vzg ;
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

  if(argc<15)
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
  sscanf(argv[12],"%lf",&M_STAR) ;
  pnmfile = argv[13];
  sscanf(argv[14],"%lu",&seed) ; // This seed is for choosing number of galaxies

  {
    unsigned long seed1 = 1234; // This seed is for choosing radial positions
    unsigned long seed2 = 9876543210; // This seed is for choosing angular positions
    work = init_nfw_workspace( seed1, seed2, NMAX_NFW );
  }
  rng = rng_create(seed);

  iarg=15 ;

  fp1=fopen(pnmfile,"w") ;

/*---Test Arguments ---*/
  if(logMmin <= 0 || logM0 <= 0 || logM1 <= 0) {
    fprintf(stderr,"halobias> Bad values for logMass!  (logMmin, logM0, logM1) = (%g,%g,%g)\n", logMmin, logM0, logM1) ;
    return(-1) ;
  }

// Scale  = Rmax(halo)/Rvir(halo)
  Scale = 1. ;

/*---Read First BGC2 file for header -----------------------------------------------------------*/

  bgc_file = argv[iarg] ;
  fp2=fopen(bgc_file,"r") ;
  assert( fp2 != NULL ) ;
  bgc_read_header(fp2,&hdr) ;
  fclose(fp2) ;
  Nhalos=hdr.ngroups_total ;
  Deltavir = hdr.overdensity ;
  mp=hdr.part_mass ; /* Rockstar sets this to full units, not 1e10 Msun/h */
  if (mp < 1e6) {
    warn("WARNING: assuming particle mass needs 1e10 factor for correct masses!");
    mp *= 1e10;
  }
  Lbox=hdr.box_size ;
  abox=1./(1+hdr.redshift) ;
  redshift = hdr.redshift;
  fprintf(stderr,"halobias> The redshift of this snapshot = %lf\n",redshift);
  fprintf(stderr,"halobias> SO halos with deltavir = %lf (wrt mean)\n",Deltavir);

/*---Read-group-info-in-BGC2-files-----------------------------------------------------*/

  /* we're allocating memory for ALL halos in BGC2 files */
  Mhalo=(float *) calloc(Nhalos,sizeof(float)) ;
  Rvir=(float *) calloc(Nhalos,sizeof(float)) ;
  xcen=(float *) calloc(Nhalos,sizeof(float)) ;
  ycen=(float *) calloc(Nhalos,sizeof(float)) ;
  zcen=(float *) calloc(Nhalos,sizeof(float)) ;
  vxcen=(float *) calloc(Nhalos,sizeof(float)) ;
  vycen=(float *) calloc(Nhalos,sizeof(float)) ;
  vzcen=(float *) calloc(Nhalos,sizeof(float)) ;

  i = j = 0 ;
  {
    int n,nfiles ;

    GROUP_DATA_RMPVMAX *gdata;
    BGC_VERBOSE = 0;

    nfiles = argc - iarg;
    if( nfiles != hdr.num_files )
    {
      warn( "WARNING: number of BGC2 files do not match expected (%d given, %d expected)\n", nfiles, hdr.num_files ) ;
    }
    for(n=0; n < nfiles; n++)
      {
        int k;
        /* Read interesting bits of BGC2 file */
        {
          FILE *fp;
          char * filename = argv[iarg + n] ;

          fp = fopen( filename, "r" );
          assert( fp != 0 );
          bgc_read_header( fp, &hdr );
          gdata = bgc_read_grouplist( fp, hdr );
          fclose( fp );
        }

        for(k=0; k < hdr.ngroups; k++)
          {
            GROUP_DATA_RMPVMAX gd = ( ( GROUP_DATA_RMPVMAX * ) gdata )[k];
            Mhalo[i] = gd.mass ;
            Rvir[i] = gd.radius ;
            xcen[i] = gd.pos[0] ;
            ycen[i] = gd.pos[1] ;
            zcen[i] = gd.pos[2] ;
            vxcen[i] = gd.vel[0] ;
            vycen[i] = gd.vel[1] ;
            vzcen[i] = gd.vel[2] ;
            if(Mhalo[i]>=Mmin)
	      {
		j++ ;
	      }
            i++ ;

          }

        free(gdata); /* clean up memory! */

      }
    if(i != Nhalos) {
      warn("WARNING: expected %d total halos, only read in %d\n", i, Nhalos);
      Nhalos = i;
    }
  }
  fprintf(stderr,"halobias> Side Length of Box = %g\n",Lbox) ;
  fprintf(stderr,"halobias> Total number of halos = %d\n",Nhalos) ;
  fprintf(stderr,"halobias> Number of halos more massive than Mmin = %d\n",j) ;

/*---Setup-galaxy-arrays--------------------------------------------------------*/

  Ngalmax = Nhalos * GALAXY_COUNT_FACTOR;
  xg =(float *) calloc(Ngalmax,sizeof(float)) ;
  yg =(float *) calloc(Ngalmax,sizeof(float)) ;
  zg =(float *) calloc(Ngalmax,sizeof(float)) ;
  vxg=(float *) calloc(Ngalmax,sizeof(float)) ;
  vyg=(float *) calloc(Ngalmax,sizeof(float)) ;
  vzg=(float *) calloc(Ngalmax,sizeof(float)) ;
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

      fprintf(fp1,"%6d %6.3f  %1d %3d\n",(halo+1),logMh,Ncen,Nsat) ;
      fflush(fp1) ;


/*---Force-a-galaxy-at-halo-center----------------------------------------------*/

      if(Ncen>0)
	{
	  if(center==1)
	    {
	      xg[gal] = (float)xcen[halo] ;
	      yg[gal] = (float)ycen[halo] ;
	      zg[gal] = (float)zcen[halo] ;
	      vxg[gal] = (float)(vxcen[halo]*sqrt(abox)) ;
	      vyg[gal] = (float)(vycen[halo]*sqrt(abox)) ;
	      vzg[gal] = (float)(vzcen[halo]*sqrt(abox)) ;
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

	      vxg[gal] = (float)(vxcen[halo]*sqrt(abox)) ;
	      vyg[gal] = (float)(vycen[halo]*sqrt(abox)) ;
	      vzg[gal] = (float)(vzcen[halo]*sqrt(abox)) ;

	      gal++;
	    }

	}
    }
  free_nfw_workspace(&work);

/*---Write-galaxy-file----------------------------------------------------------*/

  Ngal = gal ;
  fprintf(stderr,"halobias> Writing galaxy file. Ngal = %d\n",Ngal) ;

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
  ftwrite(xg,sizeof(float),Ngal,stdout);
  ftwrite(yg,sizeof(float),Ngal,stdout);
  ftwrite(zg,sizeof(float),Ngal,stdout);
  ftwrite(vxg,sizeof(float),Ngal,stdout);
  ftwrite(vyg,sizeof(float),Ngal,stdout);
  ftwrite(vzg,sizeof(float),Ngal,stdout);

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
