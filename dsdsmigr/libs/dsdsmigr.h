#ifndef __DSDSMIGR_H
#define __DSDSMIGR_H

int soho_ephemeris (TIME obs_time, double *r, double *b, double *l,
                    double *vr, double *vn, double *vw, TIME *table_mod_time);

#endif
