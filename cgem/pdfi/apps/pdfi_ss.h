/*
This include file contains function prototypes for all the pdfi_ss subroutines,
declared as void functions for calling from a C program.
c
c - -  PDFI_SS electric field inversion software
c - -  http://cgem.ssl.berkeley.edu/~fisher/public/software/PDFI_SS
c - -  Copyright (C) 2015-2019 Regents of the University of California
c 
c - -  This software is based on the concepts described in Kazachenko et al.
c - -  (2014, ApJ 795, 17).  A detailed description of the software is in
c - -  Fisher et al. (2019, arXiv:1912.08301 ).
c - -  If you use the software in a scientific 
c - -  publication, the authors would appreciate a citation to these papers 
c - -  and any future papers describing updates to the methods.
c
c - -  This is free software; you can redistribute it and/or
c - -  modify it under the terms of the GNU Lesser General Public
c - -  License as published by the Free Software Foundation,
c - -  version 2.1 of the License.
c
c - -  This software is distributed in the hope that it will be useful,
c - -  but WITHOUT ANY WARRANTY; without even the implied warranty of
c - -  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
c - -  See the GNU Lesser General Public License for more details.
c
c - -  To view the GNU Lesser General Public License visit
c - -  http://www.gnu.org/copyleft/lesser.html
c - -  or write to the Free Software Foundation, Inc.,
c - -  59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
c
*/
void abcd2wcs_ss_(int *m,int *n, double *a,double *b,double *c,double *d,
double *crpix1,double *crpix2,double *crval1,double *crval2, double *cdelt1,
double *cdelt2);
void add_padding_as_ss_(int *m,int *n, int *mpadb,int *mpadt,int *npadl,
int *npadr, double *arr, double *arrpadded);
void add_padding_ss_(int *m, int *n, int *npadlat, int *npadlon, double *arr, 
double* arrpadded);
void ahpot_ss_(int *m, int* n, int* p,int *bcn, double *a, double *b, 
double *c, double *d, double *rsun, double *rssmrs, double *scrb3d, double 
*mflux, double *atpot, double *appot);
void ahpottp2ll_ss_(int *m,int *n,int *p, double *atpepot, double *aptepot,
double *alonpot, double *alatpot);
void angle_be_ss_(int *m, int *n, double *et, double *ep, double *er, 
double *btcoe, double *bpcoe, double *brcoe, double *maskcoe, double *cosang);
void berciktest_ss_(int *m, int *n, int *p, int *bcn, double *a, double *b,
double *c, double *d, double *rsun, double *rssmrs, double *scrb3d, 
double *brh, double *brr);
void bhll2tp_ss_(int *m, int *n, double *bloncoe, double *blatcoe, 
double *btcoe, double *bpcoe);
void bhpot_phot_ss_(int *m, int *n, int *p, int *bcn, double *a, double *b,
double *c, double *d, double *rsun, double *rssmrs, double *brce, 
double *scrb3d, double *btphot, double *bpphot);
void bhpot_ss_(int *m,int *n, int *p, int *bcn, double *a,double *b, double *c,
double *d, double *rsun, double *rssmrs, double *scrb3d, double *btpot,
double *bppot);
void bhpottp2ll_ss_(int *m, int *n, int *p, double *btpot, double *bppot,
double *blonpot, double *blatpot);
void bhtp2ll_ss_(int *m,int *n, double *btcoe, double *bpcoe, double *bloncoe,
double *blatcoe);
void bhyeell2tp_ss_(int *m, int *n, double *blon, double *blat,
double *btte, double *bppe);
void bhyeetp2ll_ss_(int *m, int *n, double *btte, double *bppe, double *blon,
double *blat);
void bpot_psi_ss_(int *m, int *n,int *p, int *bcn, double *a, double*b,
double *c, double *d,double *rsun,double *rssmrs, double *psi3d, double *mflux,
double *btpot, double *bppot, double *brpot);
void br_voxels3d_ss_(int *m, int *n, double *rsun, double *sinth,
double *sinth_hlf, double *dtheta, double *dphi, double *bt, double *bp,
double *br, double *dr, double *brtop, double *brbot);
void brll2tp_ss_(int *m, int *n, double *brllcoe, double *brtpcoe);
void brpot_ss_(int *m, int *n, int *p, int *bcn, double *a, double *b,
double *c, double *d, double *rsun, double *rssmrs, double *scrb3d,
double *mflux, double *brpot);
void brpottp2ll_ss_(int *m,int *n, int *p, double *brpottp, double *brpotll);
void brtp2ll_ss_(int *m, int *n, double *brtpcoe, double *brllcoe);
void bryeell2tp_ss_(int *m, int *n, double *brllce, double *brtpce);
void bryeetp2ll_ss_(int *m, int *n, double *brtpce, double *brllce);
void bspline_ss_(double *fdata, int *nx, int *ny, double *xnew, int *nxinterp, 
double *ynew, int *nyinterp, double *finterp, int *degree);
void car2sph_ss_(int *m, int *n, double *delth, double *delx, double *dely,
double *rsph, double *a, double *b, double *c, double *d);
void cell_ss_(int *m, int *n, double *a, double *b, double *c, double *d,
double *lon, double *lat);
void cetp_ss_(int *m, int *n, double *a,double *b, double *c, double *d,
double *theta, double *phi);
void coell_ss_(int *m,int *n, double *a, double *b, double *c, double *d,
double *lon, double *lat);
void coetp_ss_(int *m,int *n, double *a, double *b, double*c, double *d,
double *theta, double *phi);
void coll_ss_(int *m, int*n, double *a, double *b, double *c, double *d,
double *lon, double *lat);
void cotp_ss_(int *m, int *n, double *a, double *b, double *c, double *d,
double *theta, double *phi);
void curl_psi_rhat_ce_ss_(int *m, int *n, double *psi, double *rsun,
double *sinth_hlf, double *dtheta, double *dphi, double *curlt, double *curlp);
void curl_psi_rhat_co_ss_(int *m, int *n, double *psi, double *rsun, 
double *sinth, double *dtheta, double *dphi, double *curlt, double *curlp);
void curlahpot_ss_(int *m, int *n, int*p, double *a, double *b, double *c,
double *d, double *rsun, double *rssmrs, double *atpot, double *appot,
double *btpot, double *bppot, double *brpot);
void curle3d_ll_(int *m, int *n, double *a, double *b, double *c, double *d,
double *rsun, double *dr, double *elontop, double *elonbot, double *elattop,
double *elatbot, double *erll, double *blont, double *blatt, double *brtll);
void curle3d_ss_(int *m,int *n, double *a, double *b, double *c, double *d,
double *rsun, double *dr, double *ettop, double *etbot, double *eptop, 
double *epbot, double *er, double *btt, double *bpt, double *brt);
void curle3dphot_ss_(int *m, int *n, double *a, double*b, double*c, double*d,
double *rsun, double *et, double *ep, double *er, double *detdr, double *depdr,
double *btt, double *bpt, double *brt);
void curle_ss_(int *m, int *n, double *rsun, double*a, double *b, double *c,
double *d, double *et, double *ep, double *scrj, double *dscrbdr,double *btt,
double *bpt, double *brt);
void curlehr_ss_(int *m, int *n, double *rsun, double *a, double *b, double *c,
double *d, double *et, double *ep, double *brt);
void curlh_ce_ss_(int *m, int *n, double *et, double *ep, double *rsun,
double *sinth, double *sinth_hlf, double *dtheta, double *dphi, double *curl);
void curlh_co_ss_(int *m, int*n, double *bt, double *bp, double *rsun,
double *sinth, double *sinth_hlf, double *dtheta, double *dphi, double *curl);
void dehdr_ss_(int *m, int*n, double*rsun, double *sinth_hlf, double *dtheta,
double *dphi, double *scrb, double *dscrbdr, double *detdr, double *depdr);
void delh2_ce_ss_(int *m, int *n, double *a, double *b, double *c, double *d,
double *rsun, double *scrb, double *delh2);
void delh2_co_ss_(int *m, int *n, double *a, double *b, double *c, double *d,
double *rsun, double *psi, double *delh2);
void delh2_sc_(int *m, int *n, double *psi, double *rsun, double *sinth_gh,
double *sinth_hlf_gh, double *dtheta, double *dphi, double *delh2);
void dilate_ss_(int *m, int *n, int *map, int *dilmap, int *dilation_param);
void divh_ce_ss_(int *m, int*n, double *bt, double *bp, double *rsun,
double *sinth, double *sinth_hlf, double*dtheta, double *dphi, double *div);
void divh_co_ss_(int *m, int *n, double *et, double *ep, double *rsun,
double *sinth, double *sinth_hlf, double *dtheta, double *dphi, double *div);
void divh_sc_(int *m, int *n, double *et, double *ep, double *rsun, 
double *sinth, double *dtheta, double *dphi, double *div);
void downsample3d_ll_(int *m, int *n, double *elon, double *elat, double *erll,
double *delondr, double *delatdr, int *mc, int*nc, double *elonc, 
double *elatc, double *erllc, double *delondrc, double *delatdrc);
void downsample3d_ss_(int *m, int*n, double *et, double *ep, double *er,
double *detdr, double *depdr, int *mc, int *nc, double *etc, double *epc,
double *erc, double *detdrc,double *depdrc);
void downsample_ll_(int *m,int *n, double *elon,double *elat,int *mc,int *nc,
double *elonc, double *elatc);
void downsample_ss_(int *m, int *n, double *et, double *ep,int *mc, int*nc,
double *etc, double *epc);
void e_doppler_rpils_ss_(int *m, int *n, double *btcoe, double *bpcoe,
double *brcoe, double *vlos, double *lt,double *lp, double *lr, double *rsun,
double *sinth, double *a, double *b, double *c, double *d, double *edopt,
double *edopp, double *edopr);
void e_doppler_ss_(int *m, int *n, double *btcoe, double *bpcoe, double *brcoe,
double *vlos,double *lt, double *lp, double*lr, double *rsun, double *sinth,
double *a, double*b, double*c, double *d, double *edopt, double *edopp,
double *edopr);
void e_flct_ss_(int *m, int*n, double *brte, double *brpe, double *bhte,
double *bhpe, double *vt, double*vp, double *rsun, double *sinth, 
double*sinth_hlf, double *a, double *b, double *c, double *d, double *ezetat,
double *ezetap);
void e_ideal_ss_(int *m, int *n, double *btcoe, double *bpcoe, double *brcoe,
double *et, double *ep, double *er, double *rsun, double *sinth, double *a,
double *b, double *c, double *d, double *eti, double *epi, double *eri);
void e_laplace_ll_(int *m,int *n, double *a, double *b, double *c, double *d,
double *rsun, double *en, double *es, double *el, double *er, double *elon,
double *elat);
void e_laplace_ss_(int *m, int*n, double *a, double *b, double *c, double *d,
double *rsun, double *ea, double *eb, double *ec, double *ed, double *et,
double *ep);
void e_ptd_ss_(int *m, int *n, double *scrbt, double *scrjt, double *rsun,
double *sinth_hlf, double *dtheta, double *dphi, double *et, double *ep,
double *er);
void e_voxels3d_ss_(int *m, int *n, double *rsun, double *sinth_hlf, 
double *dtheta, double *dphi, double *et, double *ep, double *scrb,
double *dscrbdr, double *dr, double *ettop, double *etbot, double *eptop,
double *epbot);
void ehyeell2tp_ss_(int *m,int *n, double *elon, double *elat,
double *etpe, double *epte);
void ehyeetp2ll_ss_(int *m, int*n, double *etpe, double *epte,
double *elon, double *elat);
void emagpot_psi_ss_(int *m,int *n, double *a, double *b, double *c,
double *d, double *rsun, double *psi, double *br, double *emag);
void emagpot_srf_ss_(int *m,int *n, double *a, double *b, double *c,
double *d, double *rsun, double *at, double *ap, double *bt, double *bp,
double *emag);
void emagpot_ss_(int *m, int*n, int*p, double *a, double *b, double *c,
double *d, double *rsun, double *rssmrs, double *btpot, double *bppot,
double *brpot, double *emag);
void enudge3d_gl_ll_(int *m, int*n, double *rsun, double *dr, double *blont,
double *blatt, double *brllt, double *elontop, double *elonbot, double *elattop,
double *elatbot, double *erll);
void enudge3d_gl_ss_(int *m, int *n, double *rsun, double *dr, double *btt,
double *bpt, double *brt, double *ettop, double *etbot, double *eptop,
double *epbot, double *er);
void enudge3d_ss_(int *m, int *n, double *a, double *b, double *c,
double *d, double *rsun, double *dr, double *btt, double *bpt, double *brt,
double *ettop, double *etbot, double *eptop, double *epbot, double *er);
void enudge_gl_ll_(int *m, int *n, double *rsun, double *brt, double *elon,
double *elat);
void enudge_gl_ss_(int *m, int *n, double *rsun, double *brt, double *et,
double *ep);
void enudge_ll_(int *m, int *n, double *a, double *b, double *c, double *d,
double *rsun, double *brtll, double *elon, double *elat);
void enudge_ss_(int *m, int *n, double *a, double *b, double *c, double *d,
double *rsun, double *brt, double *et, double *ep);
void eryeell2tp_ss_(int *m, int *n, double *erllcoe, double *ertpcoe);
void eryeetp2ll_ss_(int *m, int *n, double *ertpcoe, double *erllcoe);
void find_mask_ss_(int *m, int *n, double *bmag0, double *bmag1, double *bmag2,
double *bthr, double *mask);
void fix_mask_ss_(int *mdim, int *ndim, int *flag, double *mask);
void fluxbal_ll_(int *m,int *n, double *a, double *b, double *c, double *d,
double *rsun, double *brtll, double *brtllbal);
void fluxbal_ss_(int *m, int *n, double *a, double *b, double *c, double *d,
double *rsun, double *brt, double *brtbal);
void get_pils_rad_ss_(int *m, int *n, double *brad, double *bmag, int *pilmap,
double *thresh_brad, double *thresh_bmag, int *dilation_param);
void get_pils_ss_(int *m, int *n, double *bmap, int *pilmap,
double *thresh, int *dilation_param);
void gradh_ce_ss_(int *m, int *n, double *psi, double *rsun, double *sinth_hlf,
double *dtheta, double *dphi, double *gradt, double *gradp);
void gradh_co_ss_(int *m, int *n, double *psi, double *rsun, double *sinth,
double *dtheta, double *dphi, double *gradt, double *gradp);
void gradh_sc_(int *m, int *n, double *psi, double *rsun, double *sinth_gh,
double *dtheta, double *dphi, double *gradt, double *gradp);
void hm_ss_(int *m, int *n, double *rsun, double *sinth_hlf, double *dtheta,
double *dphi, double *et, double *ep, double *scrb, double *hm);
void hmtot_ss_(int *m, int *n, double *rsun, double *sinth_hlf, double *dtheta,
double *dphi, double *hm, double *hmtot);
void interp_data_ss_(int *m, int *n, double *btcoe, double *bpcoe, 
double *brcoe, double *vtcoe, double *vpcoe, double *a, double *b, double *bt,
double *bp, double *br, double *brte, double *brpe, double *bhte, double *bhpe,
double *vt, double *vp);
void interp_eh_ss_(int *m, int *n, double *et, double *ep, double *etc,
double *epc);
void interp_ehcoe_ss_(int *m, int *n, double *etcoe, double *epcoe, double *et,
double *ep);
void interp_hmidata_3d_ll_(int *m, int* n, double *data3d, int *m_int,
int *n_int, double *data3d_int);
void interp_hmidata_ll_(int *m,int *n, double *data2d, int *m_new, int *n_new,
double *data2d_new, int *degree);
void interp_var_ss_(int *m, int *n, double *btcoe, double *bpcoe, 
double *brcoe, double *a, double *b, double *bt, double *bp, double *br);
void kcost_ss_(int *n, double *k);
void kfft_ss_(int *n, double *k);
void laplacetest_ss_(int *m, int *n, int *p, int *bcn, double *a, double *b,
double *c, double *d, double *rsun, double *rssmrs, double *psi3d, double *brh,
double *brr);
void mflux_ll_(int *m, int *n, double *a, double *b,
double *c, double *d, double *rsun, double *brll, double *mflux);
void mflux_ss_(int *m,int *n, double *a, double *b, double *c, double *d,
double *rsun, double *brce, double *mflux);
void pad_abcd_as_ss_(int *m, int *n, int *mpadb, int *mpadt,int *npadl,
int *npadr, double *a, double *b, double *c, double *d, int *mp, int *np,
double *ap, double *bp, double *cp, double *dp);
void pad_abcd_ss_(int *m,int *n, int *npadlat, int *npadlon, double *a,
double *b, double *c, double *d, int *mp, int *np, double *ap, double *bp,
double *cp, double *dp);
void pad_int_gen_ss_(int *m_orig, int *n_orig, int *mpad0, int *npad0,
int *mpadb, int *mpadt, int *npadl, int *npadr);
void pdfi_wrapper4anmhd_ss_(int *m, int *n, double *rsun, double *a,double *b,
double *c, double *d, double *bloncoe0, double *blatcoe0, double *brllcoe0,
double *bloncoe1, double *blatcoe1, double *brllcoe1, double *vloncoe0,
double *vlatcoe0, double *vlosllcoe0, double *vloncoe1, double *vlatcoe1, 
double *vlosllcoe1, double *lloncoe0, double *llatcoe0, double *lrllcoe0,
double *lloncoe1, double *llatcoe1, double *lrllcoe1, double *tjul0,
double *tjul1, double *blon0, double *blat0, double *brll0, double *blon1,
double *blat1, double *brll1, double *elonpdfi, double *elatpdfi, 
double *erllpdfi, double *srll, double *srtot, double *hmll, double *hmtot,
double *tjulhalf);
void pdfi_wrapper4jsoc_ss_(int *m, int *n, double *rsun, double *a,
double *b, double *c, double *d, double *bmin, double *bloncoe0,
double *blatcoe0, double *brllcoe0, double *bloncoe1, double *blatcoe1,
double *brllcoe1, double *vloncoe0, double *vlatcoe0, double *vlosllcoe0,
double *vloncoe1, double *vlatcoe1, double *vlosllcoe1, double *lloncoe0,
double *llatcoe0, double *lrllcoe0, double *lloncoe1, double *llatcoe1,
double *lrllcoe1, double *tjul0, double *tjul1, double *blon0, double *blat0,
double *brll0, double *blon1, double *blat1, double *brll1, double *elonpdfi,
double *elatpdfi, double *erllpdfi, double *erllind, double *delondr,
double *delatdr, double *srll, double *srtot, double *hmll, double *hmtot,
double *mcoell, double *mcoll, double *mcell, double *mtell, double *mpell,
double *tjulhalf);
void pell_ss_(int *m, int *n, double *a, double *b, double *c, double *d,
double *lon, double *lat);
void petp_ss_(int *m,int *n, double *a, double *b, double *c, double *d,
double *theta, double *phi);
void psi2scrb_ss_(int *m,int *n, int *p, double *rsun, double *rssmrs,
double *psi3d, double *scrb3d);
void psi_fix_ss_(int *m, int *n, int *p, int *bcn, double *a, double *b,
double *c, double *d, double *rsun, double *rssmrs, double *psi3d,
double *mflux);
void psi_gradphot_ss_(int *m, int *n, int *p, int *bcn, double *a, double *b,
double *c, double *d, double *rsun, double *psi3d, double *btphot,
double *bpphot, double *divbh);
void psipot_phot_ss_(int *m, int *n, int *p, int *bcn, double *a, double *b,
double *c, double *d, double *rsun, double *rssmrs, double *brce, 
double *scrb3d, double *psiphot);
void psipot_ss_(int *m, int *n, int *p, int *bcn, double *a, double *b,
double *c, double *d, double *rsun, double *rssmrs, double *btte,
double *bppe, double *psi3d);
void ptdsolve_eb0_ss_(int *m, int *n, double *bt, double *bp, double *br,
double *rsun, double *sinth, double *sinth_hlf, double *a, double *b,
double *c, double *d, double *scrb, double *dscrbdr, double *scrj);
void ptdsolve_ss_(int *m, int *n, double *bt, double *bp, double *br, 
double *rsun, double *sinth, double *sinth_hlf, double *a,
double *b, double *c, double *d, double *scrb, double *dscrbdr, double *scrj);
void relax_psi_3d_ss_(int *m,int *n, double *bvec, double *evec, double *rsun,
double *a, double *b, double *c, double *d, int *max_iter, double *bthr, 
int *verbose,double *psi, double *dpsi_dr, double *etot);
void scrbpot_ss_(int *m,int *n, int *p, int *bcn, double *a, double *b, 
double *c, double *d, double *rsun, double *rssmrs, double *brce,
double *scrb3d);
void sinthta_sc_(double *thmin, double *thmax, int *m, double *sinth,
double *sinth_hlf, double *sinth_gh, double *sinth_hlf_gh);
void sinthta_ss_(double *thmin, double *thmax, int *m, double *sinth,
double *sinth_hlf);
void sr_ss_(int *m, int *n, double *et, double *ep, double *bt, double *bp,
double *sr);
void srtot_ss_(int *m, int *n, double *rsun, double *sinth_hlf, double *dtheta,
double *dphi, double *sr, double *srtot);
void stack3d_ll_(int *m, int *n, double *bloncoe0,double *blatcoe0,
double *brllcoe0, double *bloncoe1, double *blatcoe1, double *brllcoe1,
double *vloncoe0, double *vlatcoe0, double *vlosllcoe0, double *vloncoe1,
double *vlatcoe1, double *vlosllcoe1, double *lloncoe0, double *llatcoe0,
double *lrllcoe0, double *lloncoe1, double *llatcoe1, double *lrllcoe1, 
double *data_ll_3d);
void tell_ss_(int *m, int *n, double *a, double *b, double *c, double *d,
double *lon, double *lat);
void tetp_ss_(int *m, int *n, double *a, double *b, double *c, double *d,
double *theta, double *phi);
void v_ideal_ss_(int *m, int *n, double *btcoe, double *bpcoe, double *brcoe,
double *et, double *ep, double *er, double *maskcoe, double *vtcoe, 
double *vpcoe, double *vrcoe);
void wcs2abcd_ss_(int *m, int *n, double *crval1, double *crval2, 
double *cdelt1, double *cdelt2, double *a, double *b, double *c, double *d);
void wcs2mn_ss_(double *crpix1, double *crpix2, int *m, int *n);
void divb3d_ss_(int *m, int *n, int *p, double *a, double *b, double *c,
double *d, double *rsun, double *rssmrs, double *btpot, double *bppot,
double *brpot, double *divb);
