#include "astro.h"
#include "jsoc_main.h"

#define DEBUG 0

#define kXGCI "X_GEO"
#define kYGCI "Y_GEO"
#define kZGCI "Z_GEO"
#define kVXGCI "VX_GEO"
#define kVYGCI "VY_GEO"
#define kVZGCI "VZ_GEO"
#define kXHCI "X_HELIO"
#define kYHCI "Y_HELIO"
#define kZHCI "Z_HELIO"
#define kVXHCI "VX_HELIO"
#define kVYHCI "VY_HELIO"
#define kVZHCI "VZ_HELIO"
#define kOBSDATE "OBS_DATE"

/* Gets J2000.0 positions and velocities from cdf file.
   Units are km and km/s */
static void get_earth_ephem(double jd, double pos[6])
{
   long i;
   int targ=3;
   int ctr=11;
   double au=0.1495978706910000e09;

   pleph_(&jd,&targ,&ctr,pos);
   for (i=0;i<3;i++) pos[i]=pos[i]*au;
   for (i=3;i<6;i++) pos[i]=pos[i]*au/86400;
}

LIBASTRO_Error_t iorbit(DRMS_Env_t *env, const char *rsquery, LinkedList_t **info)
{
   LIBASTRO_Error_t err = kLIBASTRO_Success;
   int drmsstat = DRMS_SUCCESS;

   if (info)
   {
      char *query = strdup(rsquery);
      DRMS_RecordSet_t *rs = drms_open_records(env, query, &drmsstat);

      if (rs && rs->n > 0)
      {
         
      }

      if (query)
      {
         free(query);
      }
   }
   else
   {
      err = kLIBASTRO_InvalidArgs;
   }
   
   return err;
}

LIBASTRO_Error_t testiorbit(DRMS_Env_t *env, const char *rsquery)
{
   LIBASTRO_Error_t err = kLIBASTRO_Success;
   int drmsstat = DRMS_SUCCESS;
   char *query = strdup(rsquery);
   DRMS_RecordSet_t *rset = drms_open_records(env, query, &drmsstat);

   if (rset && rset->n > 0)
   {
      /* Get GCI values */
      DRMS_Record_t *rec = NULL;
      int irec;
      TIME obsdate;
      double **gcipos = malloc(sizeof(double *) * rset->n);
      double **gcivel = malloc(sizeof(double *) * rset->n);
      double **hcipos = malloc(sizeof(double *) * rset->n);
      double **hcivel = malloc(sizeof(double *) * rset->n);

      double au=0.1495978706910000e09; /* From JPL ephemeris */
      double pi=3.14159265358979e0;
      double deg2rad=pi/180;
      double alpha=16.13e0*deg2rad;
      double delta=26.13e0*deg2rad;
      double car0=1650.0e0;
      double carrate=4.2434255e-7; /* rotations/sec */

      /* Values for ecliptic coords */
      alpha = 75.76e0*deg2rad;
      delta = 7.25e0*deg2rad;
    
      double gcitdt,jd,help[6];
      double ex,ey,ez,evx,evy,evz; /* earth relative to sun */
      double sx,sy,sz,svx,svy,svz; /* s/c relative to sun */
      double txx,txy,txz,tyx,tyy,tyz,tzx,tzy,tzz;
      double lapp,w,rx,cx,lx,bx,vx,vy,vz;
      long icar;

      double *bs = malloc(sizeof(double) * rset->n);
      double *ls = malloc(sizeof(double) * rset->n);
      double *cs = malloc(sizeof(double) * rset->n);
      double *rs = malloc(sizeof(double) * rset->n);

      double coordrot = 0; /* angular diff between equatorial and ecliptic */
      double ey2;
      double ez2;
      double evy2;
      double evz2;

      txx=cos(alpha);
      txy=sin(alpha);
      txz=0.0;
      tyx=-sin(alpha)*cos(delta);
      tyy=cos(alpha)*cos(delta);
      tyz=sin(delta);
      tzx=sin(alpha)*sin(delta);
      tzy=-cos(alpha)*sin(delta);
      tzz=cos(delta);

      for (irec = 0; irec < rset->n; irec++)
      {
         rec = rset->records[irec];

         gcipos[irec] = malloc(sizeof(double) * 3);
         gcivel[irec] = malloc(sizeof(double) * 3);
         hcipos[irec] = malloc(sizeof(double) * 3);
         hcivel[irec] = malloc(sizeof(double) * 3);

         obsdate = drms_getkey_time(rec, kOBSDATE, NULL);

         gcipos[irec][0] = drms_getkey_double(rec, kXGCI, NULL);
         gcipos[irec][1] = drms_getkey_double(rec, kYGCI, NULL);
         gcipos[irec][2] = drms_getkey_double(rec, kZGCI, NULL);
         gcivel[irec][0] = drms_getkey_double(rec, kVXGCI, NULL);
         gcivel[irec][1] = drms_getkey_double(rec, kVYGCI, NULL);
         gcivel[irec][2] = drms_getkey_double(rec, kVZGCI, NULL);

         hcipos[irec][0] = drms_getkey_double(rec, kXHCI, NULL);
         hcipos[irec][1] = drms_getkey_double(rec, kYHCI, NULL);
         hcipos[irec][2] = drms_getkey_double(rec, kZHCI, NULL);
         hcivel[irec][0] = drms_getkey_double(rec, kVXHCI, NULL);
         hcivel[irec][1] = drms_getkey_double(rec, kVYHCI, NULL);
         hcivel[irec][2] = drms_getkey_double(rec, kVZHCI, NULL);

         gcitdt = obsdate + 32.184;
         jd=gcitdt/86400+2443144.5;

         get_earth_ephem(jd,help);
         /* results for the earth are for the same time as for SOHO, not
            after the appropriate time delay */
         ex=help[0];
         ey=help[1];
         ez=help[2];
         evx=help[3];
         evy=help[4];
         evz=help[5];

         /* eXX and eXXX (not x vals though) are all equatorial coordinates -
            convert to ecliptic coords. */
         coordrot = atan(evz/evy) - 23.45 * deg2rad;
         ey2 = ey * ey;
         ez2 = ez * ez;
         evy2 = evy * evy;
         evz2 = evz * evz;

         ey = sqrt(ey2 + ez2) * cos(coordrot);
         ez = sqrt(ey2 + ez2) * sin(coordrot);
         evy = sqrt(evy2 + evz2) * cos(coordrot);
         evz = sqrt(evy2 + evz2) * sin(coordrot);

         /* ecliptic coords */
         sx=ex+gcipos[irec][0];
         sy=ey+gcipos[irec][1];
         sz=ez+gcipos[irec][2];
         svx=evx+gcivel[irec][0];
         svy=evy+gcivel[irec][1];
         svz=evz+gcivel[irec][2];
         lapp=car0+carrate*gcitdt;
         w=(84.10+14.1844*(jd-2451545.0))*deg2rad;

#if DEBUG
         fprintf(stdout, "%-15s%-15s%-15s%-15s%-15s%-15s\n", 
                 "calc_evx", "eph_evx", "calc_evy", "eph_evy", "calc_evz", "eph_evz");
         fprintf(stdout, "%-15.8f%-15.8f%-15.8f%-15.8f%-15.8f%-15.8f\n", 
                 hcivel[irec][0] - gcivel[irec][0], 
                 evx, 
                 hcivel[irec][1] - gcivel[irec][1],
                 evy, 
                 hcivel[irec][2] - gcivel[irec][2],
                 evz);
         fprintf(stdout, "\n");
         fprintf(stdout, "%-20s%-20s%-20s%-20s%-20s%-20s\n", 
                 "calc_ex", "eph_ex", "calc_ey", "eph_ey", "calc_ez", "eph_ez");
         fprintf(stdout, "%-20.8f%-20.8f%-20.8f%-20.8f%-20.8f%-20.8f\n", 
                 hcipos[irec][0] - gcipos[irec][0], 
                 ex, 
                 hcipos[irec][1] - gcipos[irec][1],
                 ey, 
                 hcipos[irec][2] - gcipos[irec][2],
                 ez);
#endif

         /* Do s/c position relative to sun only. */
         rx=sqrt(sx*sx+sy*sy+sz*sz);
         rs[irec]=rx/au;
         bx=asin((tzx*sx+tzy*sy+tzz*sz)/rx);
         bs[irec]=bx;
         lx=atan2((tyx*sx+tyy*sy+tyz*sz),(txx*sx+txy*sy+txz*sz));
         ls[irec]=lx;
         cx=1.-fmod(lx-w,2*pi)/2/pi;
         icar=lapp-cx+0.5; /* round off lapp-le */
         cs[irec]=cx+icar;

         /* soho_orbit doesn't calculate solar x, y, and z - just velocities */
         vx=txx*svx+txy*svy+txz*svz;
         vy=tyx*svx+tyy*svy+tyz*svz;
         vz=tzx*svx+tzy*svy+tzz*svz;

         /* print out hci values (calculated and from MOC) */
         fprintf(stdout, "Position of s/c relative to sun (J2000.0) == calc_XX - derived from gci, fds_XX - from FDS data products\n");
         fprintf(stdout, "%-20s%-20s%-20s%-20s%-20s%-20s\n", 
                 "calc_x", "fds_x", "calc_y", "fds_y", "calc_z", "fds_z");
         fprintf(stdout, "%-20.8f%-20.8f%-20.8f%-20.8f%-20.8f%-20.8f\n", 
                 sx, hcipos[irec][0], hcipos[irec][1], sy, sz, hcipos[irec][2]);
         
         fprintf(stdout, "\n");

         fprintf(stdout, "Velocity of s/c relative to sun (J2000.0) == calc_XX - derived from gci, fds_XX - from FDS data products\n");
         fprintf(stdout, "%-20s%-20s%-20s%-20s%-20s%-20s\n", 
                 "calc_vx", "fds_vx", "calc_vy", "fds_vy", "calc_vz", "fds_vz");
         fprintf(stdout, "%-20.8f%-20.8f%-20.8f%-20.8f%-20.8f%-20.8f\n", 
                 svx, hcivel[irec][0], svy, hcivel[irec][1], svz, hcivel[irec][2]);
      }

      if (gcipos)
      {
         for (irec = 0; irec < rset->n; irec++)
         {
            if (gcipos[irec])
            {
               free(gcipos[irec]);
            }
         }

         free(gcipos);
      }
      
      if (gcivel)
      {
         for (irec = 0; irec < rset->n; irec++)
         {
            if (gcivel[irec])
            {
               free(gcivel[irec]);
            }
         }

         free(gcivel);
      }

      if (hcipos)
      {
         for (irec = 0; irec < rset->n; irec++)
         {
            if (hcipos[irec])
            {
               free(hcipos[irec]);
            }
         }

         free(hcipos);
      }

      if (hcivel)
      {
         for (irec = 0; irec < rset->n; irec++)
         {
            if (hcivel[irec])
            {
               free(hcivel[irec]);
            }
         }

         free(hcivel);
      }
   }

   if (query)
   {
      free(query);
   }

   return err;
}

