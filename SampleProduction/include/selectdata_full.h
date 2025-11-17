//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Jan 13 10:00:37 2024 by ROOT version 6.20/08
// from TTree amstreea/amstreea
// found on file: /eos/ams/group/mit/amsd69n_Li6MCB1306l1_3_6000605N0/1686572662_5.root
// *  History:
// *    20241029 - modified by ZX.Yan
//////////////////////////////////////////////////////////

#ifndef selectdata_h
#define selectdata_h
#define ROOT_STD_STRING_VIEW 1  // 告诉 ROOT 使用它自己的 string_view 实现

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "IsotopeAnalyzer.h"
#include <iostream>

// Header file for the classes stored in the TTree if any.

class selectdata {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   Bool_t          isbadrun;
   Bool_t          isreal;
   UInt_t          event;
   Int_t           version;
   Int_t           amsdn;
   UInt_t          errorb;
   UInt_t          time[2];
   Float_t         thetas;
   Float_t         phis;
   Float_t         rads;
   Float_t         yaw;
   Float_t         pitch;
   Float_t         roll;
   Int_t           nlevel1;
   Int_t           ndaqev;
   Int_t           ndaqerr;
   UInt_t          daqlen;
   Int_t           daqtyerr;
   UInt_t          daqjrmerr;
   UInt_t          daqjlen[24];
   UInt_t          daqsdlen;
   UInt_t          daqsdlenu;
   Float_t         livetime;
   Int_t           physbpatt;
   Int_t           jmembpatt;
   Int_t           physbpatt1;
   Int_t           jmembpatt1;
   Int_t           physbpatt2;
   Int_t           jmembpatt2;
   Int_t           ecalflag;
   Int_t           tofflag[2];
   Int_t           antipatt;
   Int_t           irti;
   Float_t         zenith;
   Float_t         glong;
   Bool_t          issaa;
   Float_t         glat;
   Float_t         thetam;
   Float_t         phim;
   Float_t         mcutoffi[4][2];
   Float_t         rtilf;
   Float_t         rtinev;
   Float_t         rtinerr;
   Float_t         rtintrig;
   Float_t         rtinpar;
   Int_t           rtigood;
   Float_t         rtinexl[2][2];
   Float_t         rticdif[2][3];
   Float_t         tk_l1q;
   Float_t         tk_l9q;
   Float_t         tk_l1qxy[2];
   Float_t         tk_l9qxy[2];
   Int_t           tk_l1qs;
   Int_t           tk_l9qs;
   Int_t           ibetahs;
   Int_t           itrdtracks;
   Int_t           tof_nhits;
   Int_t           tof_hsumhus;
   Float_t         betahs;
   Float_t         tof_chiscs;
   Float_t         tof_chists;
   Float_t         tof_qls[4];
   Float_t         trd_thetas;
   Float_t         trd_phis;
   Float_t         trd_coos[3];
   Float_t         trd_rqs;
   Float_t         ecal_pos[3];
   Float_t         ecal_ens;
   Float_t         ecal_dis;
   Float_t         tk_l1mds;
   Float_t         tk_l9mds;
   Float_t         tk_pos1s[9][3];
   Float_t         tk_dir1s[3];
   Float_t         tk_l1qvs;
   Float_t         tk_l9qvs;
   Float_t         tk_exqvn[2][2];
   Float_t         tk_l1qvr;
   Float_t         tk_l9qvr;
   Int_t           nparticle;
   Int_t           nbetah;
   Int_t           nbeta;
   Int_t           ntrack;
   Int_t           ntrrawcl;
   Int_t           ntrcl;
   Int_t           ntrhit;
   Int_t           ntrdtrack;
   Int_t           nrich;
   Int_t           nrichb;
   Int_t           nrichh;
   Int_t           nshow;
   Int_t           necalhit;
   Int_t           ntofclh;
   Int_t           ntrdseg;
   Int_t           ntrdrawh;
   Int_t           ntrdcl;
   Int_t           btstat_new;
   Int_t           ibeta;
   Int_t           ibetah;
   Int_t           icharge;
   Int_t           itrtrack;
   Int_t           itrdtrack;
   Int_t           ishow;
   Int_t           irich;
   Int_t           irichb;
   Int_t           iparindex;
   Float_t         show_dis;
   Float_t         ec_pos[3];
   Float_t         ec_dir[3];
   Float_t         cutoffpi[2];
   Int_t           betah2p;
   Float_t         betah2q;
   Float_t         betah2r;
   Int_t           betah2hb[2];
   Float_t         betah2ch;
   Int_t           tk_rtype;
   Int_t           tk_hitb[2];
   Int_t           tk_rz[2][2];
   Float_t         tk_rpz[2][2];
   Float_t         tk_qr;
   Float_t         tk_q[2];
   Float_t         tk_qrms[2];
   Int_t           tk_qhit[2];
   Float_t         tk_ql[9];
   Float_t         tk_ql2[9][2];
   Int_t           tk_qls[9];
   Int_t           tk_qlsn;
   Float_t         tk_qin[2][3];
   Float_t         tk_qrmn[2][3];
   Float_t         tk_qln[2][9][3];
   Int_t           tk_lid[9];
   Int_t           tk_cad[9];
   Float_t         tk_cof[9][2];
   Float_t         tk_iso[9][2];
   Float_t         tk_exql[2][3];
   Float_t         tk_exqln[2][2][3];
   Int_t           tk_exqls[2];
   Float_t         tk_exdis[2];
   Int_t           tk_exlid[2];
   Int_t           tk_excad[2];
   Float_t         tk_res[9][2];
   Float_t         tk_res1[9][2];
   Float_t         tk_fzm[2];
   Float_t         tk_rigidity1[3][3][7];
   Float_t         tk_rigidityk[4][7];
   Float_t         tk_erigidity1[3][7];
   Float_t         tk_chis1[3][3][7][3];
   Float_t         tk_dir0[3][3][4];
   Float_t         tk_rigidityi[9];
   Float_t         tk_chisi[9];
   Float_t         tk_bcor;
   Float_t         tk_cdif[2][3];
   Float_t         tk_pos[9][3];
   Float_t         tk_dir[9][3];
   Float_t         tk_hitc[9][3];
   Float_t         tk_oq[2];
   Float_t         tk_oel[9];
   Int_t           tk_ohitl[9];
   Int_t           tof_btype;
   Int_t           tof_bpatt;
   Float_t         tof_beta;
   Float_t         tof_betah;
   Float_t         tof_ebetah;
   Float_t         tof_ebetah_n;
   Int_t           tof_hsumh;
   Int_t           tof_hsumhu;
   Int_t           tof_nclhl[4];
   Int_t           tof_barid[4];
   Float_t         tof_pos[4][3];
   Int_t           tof_trapezoidedge[2];
   Int_t           tof_edge[4];
   Int_t           tof_goodgeo[2];
   Int_t           tof_pass;
   Float_t         tof_oq[4][2];
   Int_t           tof_ob[4][2];
   Int_t           tof_oncl;
   Float_t         tof_chist;
   Float_t         tof_chisc;
   Float_t         tof_chist_n;
   Float_t         tof_chisc_n;
   Int_t           tof_z;
   Float_t         tof_pz;
   Int_t           tof_zr;
   Float_t         tof_pzr;
   Int_t           tof_zud[2];
   Float_t         tof_pud[2];
   Float_t         tof_q[2];
   Float_t         tof_qrms[2];
   Float_t         tof_ql[4];
   Float_t         tof_ql1[4];
   Float_t         tof_qlr[3][4];
   Float_t         tof_tl[4];
   Float_t         tof_cres[4];
   UInt_t          tof_sbit[4][2];
   Int_t           tof_qs;
   Int_t           nitrdseg;
   Float_t         distrd[3];
   Int_t           trd_pass;
   Int_t           trd_statk;
   Int_t           trd_nhitk;
   Float_t         trd_lik[3];
   Float_t         trd_like[3];
   Int_t           trd_onhitk;
   Float_t         trd_oampk;
   Float_t         trd_qk[5];
   Int_t           trd_qnhk[3];
   Float_t         trd_ipch;
   Float_t         trd_qrmsk[3];
   Float_t         trd_rq;
   Int_t           trd_rz;
   Float_t         trd_rpz;
   Float_t         trd_amplk[20];
   Float_t         trd_pathlk[20];
   Float_t         trd_ampsk;
   Float_t         trd_pathsk;
   Int_t           trd_amphitk;
   Int_t           rich_itrtrack;
   Float_t         rich_beta[3];
   Float_t         rich_ebeta;
   Float_t         rich_ebeta1;
   Float_t         rich_pb;
   Float_t         rich_udis;
   Int_t           rich_hit;
   Int_t           rich_used;
   Int_t           rich_usedm;
   Int_t           rich_stat;
   Float_t         rich_q[2];
   Float_t         rich_width;
   Float_t         rich_npe[3];
   Bool_t          rich_good;
   Bool_t          rich_clean;
   Bool_t          rich_NaF;
   Int_t           rich_pmt;
   Int_t           rich_pmt1;
   Float_t         rich_PMTChargeConsistency;
   Float_t         rich_BetaConsistency;
   Float_t         rich_pos[3];
   Float_t         rich_pos_tkItr[3];
   Bool_t          rich_goodgeo;
   Float_t         rich_theta;
   Float_t         rich_theta_tkItr;
   Float_t         rich_phi;
   Float_t         rich_phi_tkItr;
   Float_t         rich_n[2];
   Int_t           rich_tile;
   Float_t         rich_distb;
   Int_t           rich_cstat;
   Float_t         rich_betap;
   Float_t         rich_likp;
   Float_t         rich_pbp;
   Int_t           rich_usedp;
   Int_t           rich_statp;
   Float_t         rich_qp;
   Float_t         rich_pzp[3];
   Int_t           rich_rz;
   Float_t         rich_rq;
   Float_t         rich_rpz;
   Int_t           anti_nhit;
   Int_t           itrtracks;
   Int_t           tk_hitbs[2];
   Float_t         tk_qis[2];
   Float_t         tk_qirmss[2];
   Float_t         tk_qlss[2][9];
   Float_t         tk_rigiditys[5];
   Float_t         tk_chiss[5][3];
   Float_t         tk_pos2s[9][3];
   Float_t         tk_dir2s[2][3];
   Int_t           tk_ibetahs;
   Float_t         tk_betahs;
   Float_t         tk_tofchiscs;
   Float_t         tk_tofchists;
   Int_t           tk_tofsumhus;
   Float_t         tk_tofqls[4];
   Float_t         iso_ll[3][3];
   Float_t         iso_da1[3][3];
   Float_t         iso_da2[3][3];
   Float_t         iso_nda1[3][3];
   Float_t         iso_nda2[3][3];
   Int_t           mpar;
   Float_t         mmom;
   Float_t         mch;
   Int_t           mtof_pass;
   Float_t         mevcoo[3];
   Float_t         mevdir[3];
   Float_t         mfcoo[3];
   Float_t         mfdir[3];
   Float_t         mpare[2];
   Int_t           mparp[2];
   Float_t         mparc[2];
   Float_t         msume;
   Float_t         msumc;
   Float_t         mtrdmom[20];
   Int_t           mtrdpar[20];
   Float_t         mtrdcoo[20][3];
   Int_t           ntrmccl;
   Float_t         mtrmom[9];
   Int_t           mtrpar[9];
   Float_t         mtrcoo[9][3];
   Float_t         mtrcoo1[9][3];
   Int_t           mtrpri[9];
   Int_t           mtrz[9];
   Int_t           ntofmccl;
   Float_t         mtofdep[4];
   Int_t           mtofpri[4];
   Int_t           mtofbar[4];
   Float_t         mevcoo1[23][3];
   Float_t         mevdir1[23][3];
   Float_t         mevmom1[23];

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_isbadrun;   //!
   TBranch        *b_isreal;   //!
   TBranch        *b_event;   //!
   TBranch        *b_version;   //!
   TBranch        *b_amsdn;   //!
   TBranch        *b_errorb;   //!
   TBranch        *b_time;   //!
   TBranch        *b_thetas;   //!
   TBranch        *b_phis;   //!
   TBranch        *b_rads;   //!
   TBranch        *b_yaw;   //!
   TBranch        *b_pitch;   //!
   TBranch        *b_roll;   //!
   TBranch        *b_nlevel1;   //!
   TBranch        *b_ndaqev;   //!
   TBranch        *b_ndaqerr;   //!
   TBranch        *b_daqlen;   //!
   TBranch        *b_daqtyerr;   //!
   TBranch        *b_daqjrmerr;   //!
   TBranch        *b_daqjlen;   //!
   TBranch        *b_daqsdlen;   //!
   TBranch        *b_daqsdlenu;   //!
   TBranch        *b_livetime;   //!
   TBranch        *b_physbpatt;   //!
   TBranch        *b_jmembpatt;   //!
   TBranch        *b_physbpatt1;   //!
   TBranch        *b_jmembpatt1;   //!
   TBranch        *b_physbpatt2;   //!
   TBranch        *b_jmembpatt2;   //!
   TBranch        *b_ecalflag;   //!
   TBranch        *b_tofflag;   //!
   TBranch        *b_antipatt;   //!
   TBranch        *b_irti;   //!
   TBranch        *b_zenith;   //!
   TBranch        *b_glong;   //!
   TBranch        *b_issaa;   //!
   TBranch        *b_glat;   //!
   TBranch        *b_thetam;   //!
   TBranch        *b_phim;   //!
   TBranch        *b_mcutoffi;   //!
   TBranch        *b_rtilf;   //!
   TBranch        *b_rtinev;   //!
   TBranch        *b_rtinerr;   //!
   TBranch        *b_rtintrig;   //!
   TBranch        *b_rtinpar;   //!
   TBranch        *b_rtigood;   //!
   TBranch        *b_rtinexl;   //!
   TBranch        *b_rticdif;   //!
   TBranch        *b_tk_l1q;   //!
   TBranch        *b_tk_l9q;   //!
   TBranch        *b_tk_l1qxy;   //!
   TBranch        *b_tk_l9qxy;   //!
   TBranch        *b_tk_l1qs;   //!
   TBranch        *b_tk_l9qs;   //!
   TBranch        *b_ibetahs;   //!
   TBranch        *b_itrdtracks;   //!
   TBranch        *b_tof_nhits;   //!
   TBranch        *b_tof_hsumhus;   //!
   TBranch        *b_betahs;   //!
   TBranch        *b_tof_chiscs;   //!
   TBranch        *b_tof_chists;   //!
   TBranch        *b_tof_qls;   //!
   TBranch        *b_trd_thetas;   //!
   TBranch        *b_trd_phis;   //!
   TBranch        *b_trd_coos;   //!
   TBranch        *b_trd_rqs;   //!
   TBranch        *b_ecal_pos;   //!
   TBranch        *b_ecal_ens;   //!
   TBranch        *b_ecal_dis;   //!
   TBranch        *b_tk_l1mds;   //!
   TBranch        *b_tk_l9mds;   //!
   TBranch        *b_tk_pos1s;   //!
   TBranch        *b_tk_dir1s;   //!
   TBranch        *b_tk_l1qvs;   //!
   TBranch        *b_tk_l9qvs;   //!
   TBranch        *b_tk_exqvn;   //!
   TBranch        *b_tk_l1qvr;   //!
   TBranch        *b_tk_l9qvr;   //!
   TBranch        *b_nparticle;   //!
   TBranch        *b_nbetah;   //!
   TBranch        *b_nbeta;   //!
   TBranch        *b_ntrack;   //!
   TBranch        *b_ntrrawcl;   //!
   TBranch        *b_ntrcl;   //!
   TBranch        *b_ntrhit;   //!
   TBranch        *b_ntrdtrack;   //!
   TBranch        *b_nrich;   //!
   TBranch        *b_nrichb;   //!
   TBranch        *b_nrichh;   //!
   TBranch        *b_nshow;   //!
   TBranch        *b_necalhit;   //!
   TBranch        *b_ntofclh;   //!
   TBranch        *b_ntrdseg;   //!
   TBranch        *b_ntrdrawh;   //!
   TBranch        *b_ntrdcl;   //!
   TBranch        *b_btstat_new;   //!
   TBranch        *b_ibeta;   //!
   TBranch        *b_ibetah;   //!
   TBranch        *b_icharge;   //!
   TBranch        *b_itrtrack;   //!
   TBranch        *b_itrdtrack;   //!
   TBranch        *b_ishow;   //!
   TBranch        *b_irich;   //!
   TBranch        *b_irichb;   //!
   TBranch        *b_iparindex;   //!
   TBranch        *b_show_dis;   //!
   TBranch        *b_ec_pos;   //!
   TBranch        *b_ec_dir;   //!
   TBranch        *b_cutoffpi;   //!
   TBranch        *b_betah2p;   //!
   TBranch        *b_betah2q;   //!
   TBranch        *b_betah2r;   //!
   TBranch        *b_betah2hb;   //!
   TBranch        *b_betah2ch;   //!
   TBranch        *b_tk_rtype;   //!
   TBranch        *b_tk_hitb;   //!
   TBranch        *b_tk_rz;   //!
   TBranch        *b_tk_rpz;   //!
   TBranch        *b_tk_qr;   //!
   TBranch        *b_tk_q;   //!
   TBranch        *b_tk_qrms;   //!
   TBranch        *b_tk_qhit;   //!
   TBranch        *b_tk_ql;   //!
   TBranch        *b_tk_ql2;   //!
   TBranch        *b_tk_qls;   //!
   TBranch        *b_tk_qlsn;   //!
   TBranch        *b_tk_qin;   //!
   TBranch        *b_tk_qrmn;   //!
   TBranch        *b_tk_qln;   //!
   TBranch        *b_tk_lid;   //!
   TBranch        *b_tk_cad;   //!
   TBranch        *b_tk_cof;   //!
   TBranch        *b_tk_iso;   //!
   TBranch        *b_tk_exql;   //!
   TBranch        *b_tk_exqln;   //!
   TBranch        *b_tk_exqls;   //!
   TBranch        *b_tk_exdis;   //!
   TBranch        *b_tk_exlid;   //!
   TBranch        *b_tk_excad;   //!
   TBranch        *b_tk_res;   //!
   TBranch        *b_tk_res1;   //!
   TBranch        *b_tk_fzm;   //!
   TBranch        *b_tk_rigidity1;   //!
   TBranch        *b_tk_rigidityk;   //!
   TBranch        *b_tk_erigidity1;   //!
   TBranch        *b_tk_chis1;   //!
   TBranch        *b_tk_dir0;   //!
   TBranch        *b_tk_rigidityi;   //!
   TBranch        *b_tk_chisi;   //!
   TBranch        *b_tk_bcor;   //!
   TBranch        *b_tk_cdif;   //!
   TBranch        *b_tk_pos;   //!
   TBranch        *b_tk_dir;   //!
   TBranch        *b_tk_hitc;   //!
   TBranch        *b_tk_oq;   //!
   TBranch        *b_tk_oel;   //!
   TBranch        *b_tk_ohitl;   //!
   TBranch        *b_tof_btype;   //!
   TBranch        *b_tof_bpatt;   //!
   TBranch        *b_tof_beta;   //!
   TBranch        *b_tof_betah;   //!
   TBranch        *b_tof_ebetah;   //!
   TBranch        *b_tof_ebetah_n;   //!
   TBranch        *b_tof_hsumh;   //!
   TBranch        *b_tof_hsumhu;   //!
   TBranch        *b_tof_nclhl;   //!
   TBranch        *b_tof_barid;   //!
   TBranch        *b_tof_pos;   //!
   TBranch        *b_tof_trapezoidedge;   //!
   TBranch        *b_tof_edge;   //!
   TBranch        *b_tof_goodgeo;   //!
   TBranch        *b_tof_pass;   //!
   TBranch        *b_tof_oq;   //!
   TBranch        *b_tof_ob;   //!
   TBranch        *b_tof_oncl;   //!
   TBranch        *b_tof_chist;   //!
   TBranch        *b_tof_chisc;   //!
   TBranch        *b_tof_chist_n;   //!
   TBranch        *b_tof_chisc_n;   //!
   TBranch        *b_tof_z;   //!
   TBranch        *b_tof_pz;   //!
   TBranch        *b_tof_zr;   //!
   TBranch        *b_tof_pzr;   //!
   TBranch        *b_tof_zud;   //!
   TBranch        *b_tof_pud;   //!
   TBranch        *b_tof_q;   //!
   TBranch        *b_tof_qrms;   //!
   TBranch        *b_tof_ql;   //!
   TBranch        *b_tof_ql1;   //!
   TBranch        *b_tof_qlr;   //!
   TBranch        *b_tof_tl;   //!
   TBranch        *b_tof_cres;   //!
   TBranch        *b_tof_sbit;   //!
   TBranch        *b_tof_qs;   //!
   TBranch        *b_nitrdseg;   //!
   TBranch        *b_distrd;   //!
   TBranch        *b_trd_pass;   //!
   TBranch        *b_trd_statk;   //!
   TBranch        *b_trd_nhitk;   //!
   TBranch        *b_trd_lik;   //!
   TBranch        *b_trd_like;   //!
   TBranch        *b_trd_onhitk;   //!
   TBranch        *b_trd_oampk;   //!
   TBranch        *b_trd_qk;   //!
   TBranch        *b_trd_qnhk;   //!
   TBranch        *b_trd_ipch;   //!
   TBranch        *b_trd_qrmsk;   //!
   TBranch        *b_trd_rq;   //!
   TBranch        *b_trd_rz;   //!
   TBranch        *b_trd_rpz;   //!
   TBranch        *b_trd_amplk;   //!
   TBranch        *b_trd_pathlk;   //!
   TBranch        *b_trd_ampsk;   //!
   TBranch        *b_trd_pathsk;   //!
   TBranch        *b_trd_amphitk;   //!
   TBranch        *b_rich_itrtrack;   //!
   TBranch        *b_rich_beta;   //!
   TBranch        *b_rich_ebeta;   //!
   TBranch        *b_rich_ebeta1;   //!
   TBranch        *b_rich_pb;   //!
   TBranch        *b_rich_udis;   //!
   TBranch        *b_rich_hit;   //!
   TBranch        *b_rich_used;   //!
   TBranch        *b_rich_usedm;   //!
   TBranch        *b_rich_stat;   //!
   TBranch        *b_rich_q;   //!
   TBranch        *b_rich_width;   //!
   TBranch        *b_rich_npe;   //!
   TBranch        *b_rich_good;   //!
   TBranch        *b_rich_clean;   //!
   TBranch        *b_rich_NaF;   //!
   TBranch        *b_rich_pmt;   //!
   TBranch        *b_rich_pmt1;   //!
   TBranch        *b_rich_PMTChargeConsistency;   //!
   TBranch        *b_rich_BetaConsistency;   //!
   TBranch        *b_rich_pos;   //!
   TBranch        *b_rich_pos_tkItr;   //!
   TBranch        *b_rich_goodgeo;   //!
   TBranch        *b_rich_theta;   //!
   TBranch        *b_rich_theta_tkItr;   //!
   TBranch        *b_rich_phi;   //!
   TBranch        *b_rich_phi_tkItr;   //!
   TBranch        *b_rich_n;   //!
   TBranch        *b_rich_tile;   //!
   TBranch        *b_rich_distb;   //!
   TBranch        *b_rich_cstat;   //!
   TBranch        *b_rich_betap;   //!
   TBranch        *b_rich_likp;   //!
   TBranch        *b_rich_pbp;   //!
   TBranch        *b_rich_usedp;   //!
   TBranch        *b_rich_statp;   //!
   TBranch        *b_rich_qp;   //!
   TBranch        *b_rich_pzp;   //!
   TBranch        *b_rich_rz;   //!
   TBranch        *b_rich_rq;   //!
   TBranch        *b_rich_rpz;   //!
   TBranch        *b_anti_nhit;   //!
   TBranch        *b_itrtracks;   //!
   TBranch        *b_tk_hitbs;   //!
   TBranch        *b_tk_qis;   //!
   TBranch        *b_tk_qirmss;   //!
   TBranch        *b_tk_qlss;   //!
   TBranch        *b_tk_rigiditys;   //!
   TBranch        *b_tk_chiss;   //!
   TBranch        *b_tk_pos2s;   //!
   TBranch        *b_tk_dir2s;   //!
   TBranch        *b_tk_ibetahs;   //!
   TBranch        *b_tk_betahs;   //!
   TBranch        *b_tk_tofchiscs;   //!
   TBranch        *b_tk_tofchists;   //!
   TBranch        *b_tk_tofsumhus;   //!
   TBranch        *b_tk_tofqls;   //!
   TBranch        *b_iso_ll;   //!
   TBranch        *b_iso_da1;   //!
   TBranch        *b_iso_da2;   //!
   TBranch        *b_iso_nda1;   //!
   TBranch        *b_iso_nda2;   //!
   TBranch        *b_mpar;   //!
   TBranch        *b_mmom;   //!
   TBranch        *b_mch;   //!
   TBranch        *b_mtof_pass;   //!
   TBranch        *b_mevcoo;   //!
   TBranch        *b_mevdir;   //!
   TBranch        *b_mfcoo;   //!
   TBranch        *b_mfdir;   //!
   TBranch        *b_mpare;   //!
   TBranch        *b_mparp;   //!
   TBranch        *b_mparc;   //!
   TBranch        *b_msume;   //!
   TBranch        *b_msumc;   //!
   TBranch        *b_mtrdmom;   //!
   TBranch        *b_mtrdpar;   //!
   TBranch        *b_mtrdcoo;   //!
   TBranch        *b_ntrmccl;   //!
   TBranch        *b_mtrmom;   //!
   TBranch        *b_mtrpar;   //!
   TBranch        *b_mtrcoo;   //!
   TBranch        *b_mtrcoo1;   //!
   TBranch        *b_mtrpri;   //!
   TBranch        *b_mtrz;   //!
   TBranch        *b_ntofmccl;   //!
   TBranch        *b_mtofdep;   //!
   TBranch        *b_mtofpri;   //!
   TBranch        *b_mtofbar;   //!
   TBranch        *b_mevcoo1;   //!
   TBranch        *b_mevdir1;   //!
   TBranch        *b_mevmom1;   //!

   void SetAnalyzer(AMS_Iso::IsotopeAnalyzer* analyzer);
   selectdata(TTree *tree=0);
   virtual ~selectdata();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

private:
   // 新增成员
   AMS_Iso::IsotopeAnalyzer* analyzer_ = nullptr;
   
};

#endif // selectdata_h

#ifdef selectdata_cxx
selectdata::selectdata(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/eos/ams/group/mit/amsd69n_Li6MCB1306l1_3_6000605N0/1686572662_5.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/eos/ams/group/mit/amsd69n_Li6MCB1306l1_3_6000605N0/1686572662_5.root");
      }
      f->GetObject("amstreea",tree);

   }
   Init(tree);
}

selectdata::~selectdata()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t selectdata::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t selectdata::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void selectdata::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("isbadrun", &isbadrun, &b_isbadrun);
   fChain->SetBranchAddress("isreal", &isreal, &b_isreal);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("version", &version, &b_version);
   fChain->SetBranchAddress("amsdn", &amsdn, &b_amsdn);
   fChain->SetBranchAddress("errorb", &errorb, &b_errorb);
   fChain->SetBranchAddress("time", time, &b_time);
   fChain->SetBranchAddress("thetas", &thetas, &b_thetas);
   fChain->SetBranchAddress("phis", &phis, &b_phis);
   fChain->SetBranchAddress("rads", &rads, &b_rads);
   fChain->SetBranchAddress("yaw", &yaw, &b_yaw);
   fChain->SetBranchAddress("pitch", &pitch, &b_pitch);
   fChain->SetBranchAddress("roll", &roll, &b_roll);
   fChain->SetBranchAddress("nlevel1", &nlevel1, &b_nlevel1);
   fChain->SetBranchAddress("ndaqev", &ndaqev, &b_ndaqev);
   fChain->SetBranchAddress("ndaqerr", &ndaqerr, &b_ndaqerr);
   fChain->SetBranchAddress("daqlen", &daqlen, &b_daqlen);
   fChain->SetBranchAddress("daqtyerr", &daqtyerr, &b_daqtyerr);
   fChain->SetBranchAddress("daqjrmerr", &daqjrmerr, &b_daqjrmerr);
   fChain->SetBranchAddress("daqjlen", daqjlen, &b_daqjlen);
   fChain->SetBranchAddress("daqsdlen", &daqsdlen, &b_daqsdlen);
   fChain->SetBranchAddress("daqsdlenu", &daqsdlenu, &b_daqsdlenu);
   fChain->SetBranchAddress("livetime", &livetime, &b_livetime);
   fChain->SetBranchAddress("physbpatt", &physbpatt, &b_physbpatt);
   fChain->SetBranchAddress("jmembpatt", &jmembpatt, &b_jmembpatt);
   fChain->SetBranchAddress("physbpatt1", &physbpatt1, &b_physbpatt1);
   fChain->SetBranchAddress("jmembpatt1", &jmembpatt1, &b_jmembpatt1);
   fChain->SetBranchAddress("physbpatt2", &physbpatt2, &b_physbpatt2);
   fChain->SetBranchAddress("jmembpatt2", &jmembpatt2, &b_jmembpatt2);
   fChain->SetBranchAddress("ecalflag", &ecalflag, &b_ecalflag);
   fChain->SetBranchAddress("tofflag", tofflag, &b_tofflag);
   fChain->SetBranchAddress("antipatt", &antipatt, &b_antipatt);
   fChain->SetBranchAddress("irti", &irti, &b_irti);
   fChain->SetBranchAddress("zenith", &zenith, &b_zenith);
   fChain->SetBranchAddress("glong", &glong, &b_glong);
   fChain->SetBranchAddress("issaa", &issaa, &b_issaa);
   fChain->SetBranchAddress("glat", &glat, &b_glat);
   fChain->SetBranchAddress("thetam", &thetam, &b_thetam);
   fChain->SetBranchAddress("phim", &phim, &b_phim);
   fChain->SetBranchAddress("mcutoffi", mcutoffi, &b_mcutoffi);
   fChain->SetBranchAddress("rtilf", &rtilf, &b_rtilf);
   fChain->SetBranchAddress("rtinev", &rtinev, &b_rtinev);
   fChain->SetBranchAddress("rtinerr", &rtinerr, &b_rtinerr);
   fChain->SetBranchAddress("rtintrig", &rtintrig, &b_rtintrig);
   fChain->SetBranchAddress("rtinpar", &rtinpar, &b_rtinpar);
   fChain->SetBranchAddress("rtigood", &rtigood, &b_rtigood);
   fChain->SetBranchAddress("rtinexl", rtinexl, &b_rtinexl);
   fChain->SetBranchAddress("rticdif", rticdif, &b_rticdif);
   fChain->SetBranchAddress("tk_l1q", &tk_l1q, &b_tk_l1q);
   fChain->SetBranchAddress("tk_l9q", &tk_l9q, &b_tk_l9q);
   fChain->SetBranchAddress("tk_l1qxy", tk_l1qxy, &b_tk_l1qxy);
   fChain->SetBranchAddress("tk_l9qxy", tk_l9qxy, &b_tk_l9qxy);
   fChain->SetBranchAddress("tk_l1qs", &tk_l1qs, &b_tk_l1qs);
   fChain->SetBranchAddress("tk_l9qs", &tk_l9qs, &b_tk_l9qs);
   fChain->SetBranchAddress("ibetahs", &ibetahs, &b_ibetahs);
   fChain->SetBranchAddress("itrdtracks", &itrdtracks, &b_itrdtracks);
   fChain->SetBranchAddress("tof_nhits", &tof_nhits, &b_tof_nhits);
   fChain->SetBranchAddress("tof_hsumhus", &tof_hsumhus, &b_tof_hsumhus);
   fChain->SetBranchAddress("betahs", &betahs, &b_betahs);
   fChain->SetBranchAddress("tof_chiscs", &tof_chiscs, &b_tof_chiscs);
   fChain->SetBranchAddress("tof_chists", &tof_chists, &b_tof_chists);
   fChain->SetBranchAddress("tof_qls", tof_qls, &b_tof_qls);
   fChain->SetBranchAddress("trd_thetas", &trd_thetas, &b_trd_thetas);
   fChain->SetBranchAddress("trd_phis", &trd_phis, &b_trd_phis);
   fChain->SetBranchAddress("trd_coos", trd_coos, &b_trd_coos);
   fChain->SetBranchAddress("trd_rqs", &trd_rqs, &b_trd_rqs);
   fChain->SetBranchAddress("ecal_pos", ecal_pos, &b_ecal_pos);
   fChain->SetBranchAddress("ecal_ens", &ecal_ens, &b_ecal_ens);
   fChain->SetBranchAddress("ecal_dis", &ecal_dis, &b_ecal_dis);
   fChain->SetBranchAddress("tk_l1mds", &tk_l1mds, &b_tk_l1mds);
   fChain->SetBranchAddress("tk_l9mds", &tk_l9mds, &b_tk_l9mds);
   fChain->SetBranchAddress("tk_pos1s", tk_pos1s, &b_tk_pos1s);
   fChain->SetBranchAddress("tk_dir1s", tk_dir1s, &b_tk_dir1s);
   fChain->SetBranchAddress("tk_l1qvs", &tk_l1qvs, &b_tk_l1qvs);
   fChain->SetBranchAddress("tk_l9qvs", &tk_l9qvs, &b_tk_l9qvs);
   fChain->SetBranchAddress("tk_exqvn", tk_exqvn, &b_tk_exqvn);
   fChain->SetBranchAddress("tk_l1qvr", &tk_l1qvr, &b_tk_l1qvr);
   fChain->SetBranchAddress("tk_l9qvr", &tk_l9qvr, &b_tk_l9qvr);
   fChain->SetBranchAddress("nparticle", &nparticle, &b_nparticle);
   fChain->SetBranchAddress("nbetah", &nbetah, &b_nbetah);
   fChain->SetBranchAddress("nbeta", &nbeta, &b_nbeta);
   fChain->SetBranchAddress("ntrack", &ntrack, &b_ntrack);
   fChain->SetBranchAddress("ntrrawcl", &ntrrawcl, &b_ntrrawcl);
   fChain->SetBranchAddress("ntrcl", &ntrcl, &b_ntrcl);
   fChain->SetBranchAddress("ntrhit", &ntrhit, &b_ntrhit);
   fChain->SetBranchAddress("ntrdtrack", &ntrdtrack, &b_ntrdtrack);
   fChain->SetBranchAddress("nrich", &nrich, &b_nrich);
   fChain->SetBranchAddress("nrichb", &nrichb, &b_nrichb);
   fChain->SetBranchAddress("nrichh", &nrichh, &b_nrichh);
   fChain->SetBranchAddress("nshow", &nshow, &b_nshow);
   fChain->SetBranchAddress("necalhit", &necalhit, &b_necalhit);
   fChain->SetBranchAddress("ntofclh", &ntofclh, &b_ntofclh);
   fChain->SetBranchAddress("ntrdseg", &ntrdseg, &b_ntrdseg);
   fChain->SetBranchAddress("ntrdrawh", &ntrdrawh, &b_ntrdrawh);
   fChain->SetBranchAddress("ntrdcl", &ntrdcl, &b_ntrdcl);
   fChain->SetBranchAddress("btstat_new", &btstat_new, &b_btstat_new);
   fChain->SetBranchAddress("ibeta", &ibeta, &b_ibeta);
   fChain->SetBranchAddress("ibetah", &ibetah, &b_ibetah);
   fChain->SetBranchAddress("icharge", &icharge, &b_icharge);
   fChain->SetBranchAddress("itrtrack", &itrtrack, &b_itrtrack);
   fChain->SetBranchAddress("itrdtrack", &itrdtrack, &b_itrdtrack);
   fChain->SetBranchAddress("ishow", &ishow, &b_ishow);
   fChain->SetBranchAddress("irich", &irich, &b_irich);
   fChain->SetBranchAddress("irichb", &irichb, &b_irichb);
   fChain->SetBranchAddress("iparindex", &iparindex, &b_iparindex);
   fChain->SetBranchAddress("show_dis", &show_dis, &b_show_dis);
   fChain->SetBranchAddress("ec_pos", ec_pos, &b_ec_pos);
   fChain->SetBranchAddress("ec_dir", ec_dir, &b_ec_dir);
   fChain->SetBranchAddress("cutoffpi", cutoffpi, &b_cutoffpi);
   fChain->SetBranchAddress("betah2p", &betah2p, &b_betah2p);
   fChain->SetBranchAddress("betah2q", &betah2q, &b_betah2q);
   fChain->SetBranchAddress("betah2r", &betah2r, &b_betah2r);
   fChain->SetBranchAddress("betah2hb", betah2hb, &b_betah2hb);
   fChain->SetBranchAddress("betah2ch", &betah2ch, &b_betah2ch);
   fChain->SetBranchAddress("tk_rtype", &tk_rtype, &b_tk_rtype);
   fChain->SetBranchAddress("tk_hitb", tk_hitb, &b_tk_hitb);
   fChain->SetBranchAddress("tk_rz", tk_rz, &b_tk_rz);
   fChain->SetBranchAddress("tk_rpz", tk_rpz, &b_tk_rpz);
   fChain->SetBranchAddress("tk_qr", &tk_qr, &b_tk_qr);
   fChain->SetBranchAddress("tk_q", tk_q, &b_tk_q);
   fChain->SetBranchAddress("tk_qrms", tk_qrms, &b_tk_qrms);
   fChain->SetBranchAddress("tk_qhit", tk_qhit, &b_tk_qhit);
   fChain->SetBranchAddress("tk_ql", tk_ql, &b_tk_ql);
   fChain->SetBranchAddress("tk_ql2", tk_ql2, &b_tk_ql2);
   fChain->SetBranchAddress("tk_qls", tk_qls, &b_tk_qls);
   fChain->SetBranchAddress("tk_qlsn", &tk_qlsn, &b_tk_qlsn);
   fChain->SetBranchAddress("tk_qin", tk_qin, &b_tk_qin);
   fChain->SetBranchAddress("tk_qrmn", tk_qrmn, &b_tk_qrmn);
   fChain->SetBranchAddress("tk_qln", tk_qln, &b_tk_qln);
   fChain->SetBranchAddress("tk_lid", tk_lid, &b_tk_lid);
   fChain->SetBranchAddress("tk_cad", tk_cad, &b_tk_cad);
   fChain->SetBranchAddress("tk_cof", tk_cof, &b_tk_cof);
   fChain->SetBranchAddress("tk_iso", tk_iso, &b_tk_iso);
   fChain->SetBranchAddress("tk_exql", tk_exql, &b_tk_exql);
   fChain->SetBranchAddress("tk_exqln", tk_exqln, &b_tk_exqln);
   fChain->SetBranchAddress("tk_exqls", tk_exqls, &b_tk_exqls);
   fChain->SetBranchAddress("tk_exdis", tk_exdis, &b_tk_exdis);
   fChain->SetBranchAddress("tk_exlid", tk_exlid, &b_tk_exlid);
   fChain->SetBranchAddress("tk_excad", tk_excad, &b_tk_excad);
   fChain->SetBranchAddress("tk_res", tk_res, &b_tk_res);
   fChain->SetBranchAddress("tk_res1", tk_res1, &b_tk_res1);
   fChain->SetBranchAddress("tk_fzm", tk_fzm, &b_tk_fzm);
   fChain->SetBranchAddress("tk_rigidity1", tk_rigidity1, &b_tk_rigidity1);
   fChain->SetBranchAddress("tk_rigidityk", tk_rigidityk, &b_tk_rigidityk);
   fChain->SetBranchAddress("tk_erigidity1", tk_erigidity1, &b_tk_erigidity1);
   fChain->SetBranchAddress("tk_chis1", tk_chis1, &b_tk_chis1);
   fChain->SetBranchAddress("tk_dir0", tk_dir0, &b_tk_dir0);
   fChain->SetBranchAddress("tk_rigidityi", tk_rigidityi, &b_tk_rigidityi);
   fChain->SetBranchAddress("tk_chisi", tk_chisi, &b_tk_chisi);
   fChain->SetBranchAddress("tk_bcor", &tk_bcor, &b_tk_bcor);
   fChain->SetBranchAddress("tk_cdif", tk_cdif, &b_tk_cdif);
   fChain->SetBranchAddress("tk_pos", tk_pos, &b_tk_pos);
   fChain->SetBranchAddress("tk_dir", tk_dir, &b_tk_dir);
   fChain->SetBranchAddress("tk_hitc", tk_hitc, &b_tk_hitc);
   fChain->SetBranchAddress("tk_oq", tk_oq, &b_tk_oq);
   fChain->SetBranchAddress("tk_oel", tk_oel, &b_tk_oel);
   fChain->SetBranchAddress("tk_ohitl", tk_ohitl, &b_tk_ohitl);
   fChain->SetBranchAddress("tof_btype", &tof_btype, &b_tof_btype);
   fChain->SetBranchAddress("tof_bpatt", &tof_bpatt, &b_tof_bpatt);
   fChain->SetBranchAddress("tof_beta", &tof_beta, &b_tof_beta);
   fChain->SetBranchAddress("tof_betah", &tof_betah, &b_tof_betah);
   fChain->SetBranchAddress("tof_ebetah", &tof_ebetah, &b_tof_ebetah);
   fChain->SetBranchAddress("tof_ebetah_n", &tof_ebetah_n, &b_tof_ebetah_n);
   fChain->SetBranchAddress("tof_hsumh", &tof_hsumh, &b_tof_hsumh);
   fChain->SetBranchAddress("tof_hsumhu", &tof_hsumhu, &b_tof_hsumhu);
   fChain->SetBranchAddress("tof_nclhl", tof_nclhl, &b_tof_nclhl);
   fChain->SetBranchAddress("tof_barid", tof_barid, &b_tof_barid);
   fChain->SetBranchAddress("tof_pos", tof_pos, &b_tof_pos);
   fChain->SetBranchAddress("tof_trapezoidedge", tof_trapezoidedge, &b_tof_trapezoidedge);
   fChain->SetBranchAddress("tof_edge", tof_edge, &b_tof_edge);
   fChain->SetBranchAddress("tof_goodgeo", tof_goodgeo, &b_tof_goodgeo);
   fChain->SetBranchAddress("tof_pass", &tof_pass, &b_tof_pass);
   fChain->SetBranchAddress("tof_oq", tof_oq, &b_tof_oq);
   fChain->SetBranchAddress("tof_ob", tof_ob, &b_tof_ob);
   fChain->SetBranchAddress("tof_oncl", &tof_oncl, &b_tof_oncl);
   fChain->SetBranchAddress("tof_chist", &tof_chist, &b_tof_chist);
   fChain->SetBranchAddress("tof_chisc", &tof_chisc, &b_tof_chisc);
   fChain->SetBranchAddress("tof_chist_n", &tof_chist_n, &b_tof_chist_n);
   fChain->SetBranchAddress("tof_chisc_n", &tof_chisc_n, &b_tof_chisc_n);
   fChain->SetBranchAddress("tof_z", &tof_z, &b_tof_z);
   fChain->SetBranchAddress("tof_pz", &tof_pz, &b_tof_pz);
   fChain->SetBranchAddress("tof_zr", &tof_zr, &b_tof_zr);
   fChain->SetBranchAddress("tof_pzr", &tof_pzr, &b_tof_pzr);
   fChain->SetBranchAddress("tof_zud", tof_zud, &b_tof_zud);
   fChain->SetBranchAddress("tof_pud", tof_pud, &b_tof_pud);
   fChain->SetBranchAddress("tof_q", tof_q, &b_tof_q);
   fChain->SetBranchAddress("tof_qrms", tof_qrms, &b_tof_qrms);
   fChain->SetBranchAddress("tof_ql", tof_ql, &b_tof_ql);
   fChain->SetBranchAddress("tof_ql1", tof_ql1, &b_tof_ql1);
   fChain->SetBranchAddress("tof_qlr", tof_qlr, &b_tof_qlr);
   fChain->SetBranchAddress("tof_tl", tof_tl, &b_tof_tl);
   fChain->SetBranchAddress("tof_cres", tof_cres, &b_tof_cres);
   fChain->SetBranchAddress("tof_sbit", tof_sbit, &b_tof_sbit);
   fChain->SetBranchAddress("tof_qs", &tof_qs, &b_tof_qs);
   fChain->SetBranchAddress("nitrdseg", &nitrdseg, &b_nitrdseg);
   fChain->SetBranchAddress("distrd", distrd, &b_distrd);
   fChain->SetBranchAddress("trd_pass", &trd_pass, &b_trd_pass);
   fChain->SetBranchAddress("trd_statk", &trd_statk, &b_trd_statk);
   fChain->SetBranchAddress("trd_nhitk", &trd_nhitk, &b_trd_nhitk);
   fChain->SetBranchAddress("trd_lik", trd_lik, &b_trd_lik);
   fChain->SetBranchAddress("trd_like", trd_like, &b_trd_like);
   fChain->SetBranchAddress("trd_onhitk", &trd_onhitk, &b_trd_onhitk);
   fChain->SetBranchAddress("trd_oampk", &trd_oampk, &b_trd_oampk);
   fChain->SetBranchAddress("trd_qk", trd_qk, &b_trd_qk);
   fChain->SetBranchAddress("trd_qnhk", trd_qnhk, &b_trd_qnhk);
   fChain->SetBranchAddress("trd_ipch", &trd_ipch, &b_trd_ipch);
   fChain->SetBranchAddress("trd_qrmsk", trd_qrmsk, &b_trd_qrmsk);
   fChain->SetBranchAddress("trd_rq", &trd_rq, &b_trd_rq);
   fChain->SetBranchAddress("trd_rz", &trd_rz, &b_trd_rz);
   fChain->SetBranchAddress("trd_rpz", &trd_rpz, &b_trd_rpz);
   fChain->SetBranchAddress("trd_amplk", trd_amplk, &b_trd_amplk);
   fChain->SetBranchAddress("trd_pathlk", trd_pathlk, &b_trd_pathlk);
   fChain->SetBranchAddress("trd_ampsk", &trd_ampsk, &b_trd_ampsk);
   fChain->SetBranchAddress("trd_pathsk", &trd_pathsk, &b_trd_pathsk);
   fChain->SetBranchAddress("trd_amphitk", &trd_amphitk, &b_trd_amphitk);
   fChain->SetBranchAddress("rich_itrtrack", &rich_itrtrack, &b_rich_itrtrack);
   fChain->SetBranchAddress("rich_beta", rich_beta, &b_rich_beta);
   fChain->SetBranchAddress("rich_ebeta", &rich_ebeta, &b_rich_ebeta);
   fChain->SetBranchAddress("rich_ebeta1", &rich_ebeta1, &b_rich_ebeta1);
   fChain->SetBranchAddress("rich_pb", &rich_pb, &b_rich_pb);
   fChain->SetBranchAddress("rich_udis", &rich_udis, &b_rich_udis);
   fChain->SetBranchAddress("rich_hit", &rich_hit, &b_rich_hit);
   fChain->SetBranchAddress("rich_used", &rich_used, &b_rich_used);
   fChain->SetBranchAddress("rich_usedm", &rich_usedm, &b_rich_usedm);
   fChain->SetBranchAddress("rich_stat", &rich_stat, &b_rich_stat);
   fChain->SetBranchAddress("rich_q", rich_q, &b_rich_q);
   fChain->SetBranchAddress("rich_width", &rich_width, &b_rich_width);
   fChain->SetBranchAddress("rich_npe", rich_npe, &b_rich_npe);
   fChain->SetBranchAddress("rich_good", &rich_good, &b_rich_good);
   fChain->SetBranchAddress("rich_clean", &rich_clean, &b_rich_clean);
   fChain->SetBranchAddress("rich_NaF", &rich_NaF, &b_rich_NaF);
   fChain->SetBranchAddress("rich_pmt", &rich_pmt, &b_rich_pmt);
   fChain->SetBranchAddress("rich_pmt1", &rich_pmt1, &b_rich_pmt1);
   fChain->SetBranchAddress("rich_PMTChargeConsistency", &rich_PMTChargeConsistency, &b_rich_PMTChargeConsistency);
   fChain->SetBranchAddress("rich_BetaConsistency", &rich_BetaConsistency, &b_rich_BetaConsistency);
   fChain->SetBranchAddress("rich_pos", rich_pos, &b_rich_pos);
   fChain->SetBranchAddress("rich_pos_tkItr", rich_pos_tkItr, &b_rich_pos_tkItr);
   fChain->SetBranchAddress("rich_goodgeo", &rich_goodgeo, &b_rich_goodgeo);
   fChain->SetBranchAddress("rich_theta", &rich_theta, &b_rich_theta);
   fChain->SetBranchAddress("rich_theta_tkItr", &rich_theta_tkItr, &b_rich_theta_tkItr);
   fChain->SetBranchAddress("rich_phi", &rich_phi, &b_rich_phi);
   fChain->SetBranchAddress("rich_phi_tkItr", &rich_phi_tkItr, &b_rich_phi_tkItr);
   fChain->SetBranchAddress("rich_n", rich_n, &b_rich_n);
   fChain->SetBranchAddress("rich_tile", &rich_tile, &b_rich_tile);
   fChain->SetBranchAddress("rich_distb", &rich_distb, &b_rich_distb);
   fChain->SetBranchAddress("rich_cstat", &rich_cstat, &b_rich_cstat);
   fChain->SetBranchAddress("rich_betap", &rich_betap, &b_rich_betap);
   fChain->SetBranchAddress("rich_likp", &rich_likp, &b_rich_likp);
   fChain->SetBranchAddress("rich_pbp", &rich_pbp, &b_rich_pbp);
   fChain->SetBranchAddress("rich_usedp", &rich_usedp, &b_rich_usedp);
   fChain->SetBranchAddress("rich_statp", &rich_statp, &b_rich_statp);
   fChain->SetBranchAddress("rich_qp", &rich_qp, &b_rich_qp);
   fChain->SetBranchAddress("rich_pzp", rich_pzp, &b_rich_pzp);
   fChain->SetBranchAddress("rich_rz", &rich_rz, &b_rich_rz);
   fChain->SetBranchAddress("rich_rq", &rich_rq, &b_rich_rq);
   fChain->SetBranchAddress("rich_rpz", &rich_rpz, &b_rich_rpz);
   fChain->SetBranchAddress("anti_nhit", &anti_nhit, &b_anti_nhit);
   fChain->SetBranchAddress("itrtracks", &itrtracks, &b_itrtracks);
   fChain->SetBranchAddress("tk_hitbs", tk_hitbs, &b_tk_hitbs);
   fChain->SetBranchAddress("tk_qis", tk_qis, &b_tk_qis);
   fChain->SetBranchAddress("tk_qirmss", tk_qirmss, &b_tk_qirmss);
   fChain->SetBranchAddress("tk_qlss", tk_qlss, &b_tk_qlss);
   fChain->SetBranchAddress("tk_rigiditys", tk_rigiditys, &b_tk_rigiditys);
   fChain->SetBranchAddress("tk_chiss", tk_chiss, &b_tk_chiss);
   fChain->SetBranchAddress("tk_pos2s", tk_pos2s, &b_tk_pos2s);
   fChain->SetBranchAddress("tk_dir2s", tk_dir2s, &b_tk_dir2s);
   fChain->SetBranchAddress("tk_ibetahs", &tk_ibetahs, &b_tk_ibetahs);
   fChain->SetBranchAddress("tk_betahs", &tk_betahs, &b_tk_betahs);
   fChain->SetBranchAddress("tk_tofchiscs", &tk_tofchiscs, &b_tk_tofchiscs);
   fChain->SetBranchAddress("tk_tofchists", &tk_tofchists, &b_tk_tofchists);
   fChain->SetBranchAddress("tk_tofsumhus", &tk_tofsumhus, &b_tk_tofsumhus);
   fChain->SetBranchAddress("tk_tofqls", tk_tofqls, &b_tk_tofqls);
   fChain->SetBranchAddress("iso_ll", iso_ll, &b_iso_ll);
   fChain->SetBranchAddress("iso_da1", iso_da1, &b_iso_da1);
   fChain->SetBranchAddress("iso_da2", iso_da2, &b_iso_da2);
   fChain->SetBranchAddress("iso_nda1", iso_nda1, &b_iso_nda1);
   fChain->SetBranchAddress("iso_nda2", iso_nda2, &b_iso_nda2);
   fChain->SetBranchAddress("mpar", &mpar, &b_mpar);
   fChain->SetBranchAddress("mmom", &mmom, &b_mmom);
   fChain->SetBranchAddress("mch", &mch, &b_mch);
   fChain->SetBranchAddress("mtof_pass", &mtof_pass, &b_mtof_pass);
   fChain->SetBranchAddress("mevcoo", mevcoo, &b_mevcoo);
   fChain->SetBranchAddress("mevdir", mevdir, &b_mevdir);
   fChain->SetBranchAddress("mfcoo", mfcoo, &b_mfcoo);
   fChain->SetBranchAddress("mfdir", mfdir, &b_mfdir);
   fChain->SetBranchAddress("mpare", mpare, &b_mpare);
   fChain->SetBranchAddress("mparp", mparp, &b_mparp);
   fChain->SetBranchAddress("mparc", mparc, &b_mparc);
   fChain->SetBranchAddress("msume", &msume, &b_msume);
   fChain->SetBranchAddress("msumc", &msumc, &b_msumc);
   fChain->SetBranchAddress("mtrdmom", mtrdmom, &b_mtrdmom);
   fChain->SetBranchAddress("mtrdpar", mtrdpar, &b_mtrdpar);
   fChain->SetBranchAddress("mtrdcoo", mtrdcoo, &b_mtrdcoo);
   fChain->SetBranchAddress("ntrmccl", &ntrmccl, &b_ntrmccl);
   fChain->SetBranchAddress("mtrmom", mtrmom, &b_mtrmom);
   fChain->SetBranchAddress("mtrpar", mtrpar, &b_mtrpar);
   fChain->SetBranchAddress("mtrcoo", mtrcoo, &b_mtrcoo);
   fChain->SetBranchAddress("mtrcoo1", mtrcoo1, &b_mtrcoo1);
   fChain->SetBranchAddress("mtrpri", mtrpri, &b_mtrpri);
   fChain->SetBranchAddress("mtrz", mtrz, &b_mtrz);
   fChain->SetBranchAddress("ntofmccl", &ntofmccl, &b_ntofmccl);
   fChain->SetBranchAddress("mtofdep", mtofdep, &b_mtofdep);
   fChain->SetBranchAddress("mtofpri", mtofpri, &b_mtofpri);
   fChain->SetBranchAddress("mtofbar", mtofbar, &b_mtofbar);
   fChain->SetBranchAddress("mevcoo1", mevcoo1, &b_mevcoo1);
   fChain->SetBranchAddress("mevdir1", mevdir1, &b_mevdir1);
   fChain->SetBranchAddress("mevmom1", mevmom1, &b_mevmom1);
   Notify();
}

Bool_t selectdata::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void selectdata::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t selectdata::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // selectdata_cxx
