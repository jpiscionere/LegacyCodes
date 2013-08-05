/* PROGRAM MAKEMOCK

   --- makemock geometry zspace ifib rotanglex rotabgley rotanglez zmin zmax maskfile < partfile > outfile
   --- construct a mock galaxy sample from a cubical galaxy distribution.
      * geometry  = 0 : sphere
                  = 1 : SDSS DR72
      * zspace    = 0 : real space
                  = 1 : redshift space 
      * ifib      = 0 : Do not impose fiber collisions
                  = 1 : Remove collided galaxies
                  = 2 : Put collided galaxies at redshift of nearest neighbors
      * rotanglex = rotation angle around x-axis (in degrees)
      * rotangley = rotation angle around y-axis (in degrees)
      * rotanglez = rotation angle around z-axis (in degrees)
      * zmin      = minimum redshift of sample
      * zmax      = maximum redshift of sample
      * maskfile  = SDSSPIX version of polygon file (pix format)
      * fgot_min  = minimum acceptable value of completeness
      * tbin_out  = 0 : ascii ra,dec,cz
      *             1 : tipsy binary format comoving x,y,z
      < partfile  = input particle file in fastfood format
      > outfile   = output file in vollim format.
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <nr.h>

#include "tbin.h"

#define sqr(x) ((x)*(x))
#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))
#define mabs(A) ((A) < 0.0 ? -(A) : (A))
#define cnint(x) ((x-floor(x)) < 0.5 ? floor(x) : ceil(x))
#define pi 3.141592653589
#define c 299800.

#define COSMO_DIST_SIZE 10000
#define OMEGA_M 0.25 
#define OMEGA_L 0.75

int main(int argc, char *argv[])
{
  int i,j,k ;
/*---Arguments-------------------------*/
  int geometry,zspace,ifib,tbin_out=0;
  double rotanglex,rotangley,rotanglez,zmin,zmax,fgot_min ;
/*---Particle-distribution-variables---*/
  int np,ibyte,*idat ;
  float *fdat,znow ;
  float Lbox ;
  float *x,*y,*z,*vx,*vy,*vz ;
/*---Cosmodistfile-variables-----------*/
  int Nzdc ;
  double *zc,*dc ;
/*---Redshift-space--------------------*/
  double r,vr,*cz1 ;
/*---SDSS-mask-------------------------*/
  int *ID,Ngal ;
  double comp1,*comp,*ra,*dec,*cz ;
  double ra1,dec1,vec[3] ;
  char *pixfile ;
/*---Fiber-collisions------------------*/
  int Nfibmax,Nfibgrp,indx,jndx,Ni,indx2,done,ipick ;
  int **IDfib,*fiblist,*idrank,*fibcollisions ;
  long seed=-1 ;
  double *x1,*y1,*z1,length ;
  double costheta,cosfiber,*costhetarank ;
/*---Functions-------------------------*/
  void InitializePolygonMask(char *) ;
  double FindCompleteness(double, double) ;
  void rotate_x(double,double *) ;
  void rotate_y(double,double *) ;
  void rotate_z(double,double *) ;
  int set_cosmo_dist(const double Omegam, const double OmegaL, const double zmax, const int maxsize, 
                     double * zc, double * dc ) ;

/*---Read-Arguments---------------------------------------------------------*/

  if(argc < 11) 
    { 
      fprintf(stderr, "Usage:\n  %s  geometry zspace ifib rotanglex rotangley rotanglez zmin zmax maskfile fgot_min tbin_out < partfile > outfile\n",
              argv[0]) ;
      exit(1) ;
    } 

  sscanf(argv[1],"%d",&geometry) ;
  sscanf(argv[2],"%d",&zspace) ;
  sscanf(argv[3],"%d",&ifib) ;
  sscanf(argv[4],"%lf",&rotanglex) ;
  sscanf(argv[5],"%lf",&rotangley) ;
  sscanf(argv[6],"%lf",&rotanglez) ;
  sscanf(argv[7],"%lf",&zmin) ;
  sscanf(argv[8],"%lf",&zmax) ;
  pixfile = argv[9];
  sscanf(argv[10],"%lf",&fgot_min) ;
  if(argc >= 12) 
    { 
      sscanf(argv[11],"%d",&tbin_out) ;
    } 

/*---Read-particle-file-----------------------------------------------------*/

  idat=(int *)calloc(5,sizeof(int));
  fdat=(float *) calloc(9,sizeof(float));       

  ftread(idat,sizeof(int),5,stdin);
  ftread(fdat,sizeof(float),9,stdin);
  ftread(&znow,sizeof(float),1,stdin) ;
 
  np=idat[1];
  Lbox=fdat[0];

  x=(float *) calloc(np,sizeof(float));       
  y=(float *) calloc(np,sizeof(float));       
  z=(float *) calloc(np,sizeof(float));       
  vx=(float *) calloc(np,sizeof(float));       
  vy=(float *) calloc(np,sizeof(float));       
  vz=(float *) calloc(np,sizeof(float));       
      
  ftread(x,sizeof(float),np,stdin);
  ftread(y,sizeof(float),np,stdin);
  ftread(z,sizeof(float),np,stdin);
  ftread(vx,sizeof(float),np,stdin);
  ftread(vy,sizeof(float),np,stdin);
  ftread(vz,sizeof(float),np,stdin);

  fprintf(stderr,"makemock> Read galaxy distribution, np=%d\n",np) ;
  fprintf(stderr,"makemock> Size of box, Lbox=%g\n",Lbox) ;

/*---Use-center-of-box-as-origin--------------------------------------------*/

  for (i=0;i<np;i++) 
    {
      x[i] -= 0.5*Lbox ;
      y[i] -= 0.5*Lbox ;
      z[i] -= 0.5*Lbox ;
    }

/*--- Calculate cosmo distance ----------------------------------------------------*/
  zc=(double *) calloc(COSMO_DIST_SIZE ,sizeof(double));
  dc=(double *) calloc(COSMO_DIST_SIZE ,sizeof(double));

  Nzdc = set_cosmo_dist(OMEGA_M, OMEGA_L, zmax, COSMO_DIST_SIZE, zc, dc);

/*---Compute-cosmological-redshifts-for-galaxies----------------------------*/
  cz1=(double *) calloc(np,sizeof(double));
  fprintf(stderr,"makemock> Computing cosmological redshifts\n") ;
  { 
      double cz_local ;
      for(i=0;i<np;i++)
        {
          cz_local = -1.0;
          r = (double)sqrt(sqr(x[i])+sqr(y[i])+sqr(z[i])) ;
          for(j=1;j<Nzdc;j++)
            {
              if(r<dc[j])
                {
                  cz_local = c*(zc[j] - (dc[j]-r)*(zc[j]-zc[j-1])/(dc[j]-dc[j-1])) ;
                  break ;
                }
            }
          if( cz_local > 0 ) { 
              cz1[i] = cz_local ;
          } else { 
              fprintf(stderr, "makemock> ERROR! couldn't set cosmological redshift for particle %d\n", i) ;
              exit(1) ;

          } 
        }
  }

/*---Put-in-redshift-space-with-observer-at-center--------------------------*/
  if(zspace==1)
    {
      for(i=0;i<np;i++)
	{
	  r = (double)sqrt(sqr(x[i])+sqr(y[i])+sqr(z[i])) ;
	  
	  vr = (double)(vx[i]*x[i]+vy[i]*y[i]+vz[i]*z[i])/r ; 

	  cz1[i] += vr + vr*cz1[i] / c; /* too many factors of c, have to remove one */
	}
      fprintf(stderr,"makemock> Put galaxies in redshift space\n") ;
    }

/*---Initialize-Mask--------------------------------------------------------*/
  if(geometry==1)
    {
      fprintf(stderr,"makemock> Trim to SDSS geometry: filtering completeness < %g\n", fgot_min) ;
      InitializePolygonMask(pixfile) ;
    }

/*---Only-keep-galaxies-that-are-in-SDSS-volume-----------------------------*/

  ID=(int *) calloc(np,sizeof(int));
  ra=(double *) calloc(np,sizeof(double));       
  dec=(double *) calloc(np,sizeof(double));       
  cz=(double *) calloc(np,sizeof(double));       
  x1=(double *) calloc(np,sizeof(double));       
  y1=(double *) calloc(np,sizeof(double));       
  z1=(double *) calloc(np,sizeof(double));       
  comp=(double *) calloc(np,sizeof(double));

  Ngal=0 ;
  for(i=0;i<np;i++)
    {
      if(cz1[i]>=zmin*c && cz1[i]<=zmax*c)
	{
	  vec[0] = (double)x[i] ;
	  vec[1] = (double)y[i] ;
	  vec[2] = (double)z[i] ;
	  length = sqrt(sqr(vec[0])+sqr(vec[1])+sqr(vec[2])) ;
	  vec[0] /= length ;
	  vec[1] /= length ;
	  vec[2] /= length ;

	  rotate_x(rotanglex,vec) ;
	  rotate_y(rotangley,vec) ;
	  rotate_z(rotanglez,vec) ;

	  dec1 = 90. - (180./pi)*acos(vec[2]) ;
	  ra1 = (180./pi)*atan(vec[1]/vec[0]) ;
	  
	  if(vec[0]<0) 
	    ra1 += 180 ;
	  else if(vec[0]>=0 && vec[1]<0) 
	    ra1 += 360 ;

	  if(geometry==0)
	    {
	      comp1=1. ;
	    }
	  else if(geometry==1)
	    {
	      comp1 = FindCompleteness(ra1,dec1) ;
	    }

	  if(comp1 >= fgot_min)
	    {
	      ID[Ngal] = Ngal ;
	      cz[Ngal] = cz1[i] ;
	      dec[Ngal] = dec1 ;
	      ra[Ngal] = ra1 ;
	      comp[Ngal] = comp1 ;

              if(zspace == 1) 
                { 
                  double dist;
                  double z_local = (cz1[i] / c);

                  for(j=1;j<Nzdc;j++)
                  {
                      if(z_local < zc[j])
                      {
                          dist =  dc[j] - (zc[j] - z_local)*(dc[j]-dc[j-1])/(zc[j]-zc[j-1]) ;
                          break ;
                      }
                  }
                  x1[Ngal] = dist*vec[0] ;
                  y1[Ngal] = dist*vec[1] ;
                  z1[Ngal] = dist*vec[2] ;
                } 
              else 
                { 
                  x1[Ngal] = (double)x[i] ; 
                  y1[Ngal] = (double)y[i] ; 
                  z1[Ngal] = (double)z[i] ; 
                }

	      Ngal++ ;
	    }
	}
    }
  fprintf(stderr,"makemock> Ngal=%d\n",Ngal) ;
  { 
    /* begin output */
    XDR xdrs;
    struct dump h;
    struct dark_particle dp;
    fprintf(stderr,"makemock> creating output in %s \n", tbin_out ? "tipsy binary format" : "ascii" );

    if(tbin_out) 
      { 
        /* Make TIPSY header */
        h.nbodies = Ngal ;
        h.ndark   = h.nbodies ;
        h.nstar   = 0 ;
        h.nsph    = 0 ;
        h.ndim    = 3 ;
        h.time    = 1.0 ;
        xdrstdio_create(&xdrs, stdout, XDR_ENCODE);
        if( xdr_header(&xdrs,&h) != 1 )
          { 
            fprintf(stderr, "Error writing tipsy binary header!");
            exit(2);
          }
        dp.mass = 1.0 ;
        for (i=0; i<3; ++i) 
          {
            dp.vel[i] = 0.0 ;
            dp.pos[i] = 0.0 ;
          }
        dp.phi = 0.0 ;
      }
    for(i=0;i<Ngal; i++) 
      { 
          if(tbin_out) 
            { 
              dp.pos[0] = x1[i] ;
              dp.pos[1] = y1[i] ;
              dp.pos[2] = z1[i] ;
              if( xdr_dark(&xdrs,&dp) != 1 )
                { 
                  fprintf(stderr, "Error writing tipsy binary data! (particle %d)",i);
                  exit(2);
                }
            } 
          else 
            { 
              fprintf(stdout,"%10.6f % 10.6f  %10.2f\n",ra[i],dec[i],cz[i]) ;
//               fprintf(stdout,"%15.6f %15.6f %15.6f \n",x1[i],y1[i],z1[i]) ;
// 	      fprintf(stdout,"%6d %10.6f %10.6f  %8.2f\n",ID[Ngal],ra[Ngal],dec[Ngal],cz[Ngal]) ;
            }
      } 
    if(tbin_out) 
      { 
        /* close XDR file */
        xdr_destroy(&xdrs);
      } 

  } 
  return 0 ;
}
/*---SUBROUTINE ROTATE_X--------------------------
--- rotate_x(double phi, double *x)
--- rotate a vector by angle phi about the X-axis
-------------------------------------------------*/
void rotate_x(double phi, double *x)
{
  double A,B ;

  A = x[1] ;
  B = x[2] ;

  x[1] = cos(phi*pi/180.)*A - sin(phi*pi/180.)*B ;
  x[2] = sin(phi*pi/180.)*A + cos(phi*pi/180.)*B ;
}
/*---SUBROUTINE ROTATE_Y--------------------------
--- rotate_y(double phi, double *x)
--- rotate a vector by angle phi about the Y-axis
-------------------------------------------------*/
void rotate_y(double phi, double *x)
{
  double A,B ;

  A = x[0] ;
  B = x[2] ;

  x[0] = cos(phi*pi/180.)*A - sin(phi*pi/180.)*B ;
  x[2] = sin(phi*pi/180.)*A + cos(phi*pi/180.)*B ;
}
/*---SUBROUTINE ROTATE_Z--------------------------
--- rotate_z(double phi, double *x)
--- rotate a vector by angle phi about the Z-axis
-------------------------------------------------*/
void rotate_z(double phi, double *x)
{
  double A,B ;

  A = x[0] ;
  B = x[1] ;

  x[0] = cos(phi*pi/180.)*A - sin(phi*pi/180.)*B ;
  x[1] = sin(phi*pi/180.)*A + cos(phi*pi/180.)*B ;
}

