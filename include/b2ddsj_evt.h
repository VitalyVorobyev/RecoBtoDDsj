#ifndef __B2DDSJ_EVT_H__
#define __B2DDSJ_EVT_H__

#include "TObject.h"
#include "mytools/datastructures.h"

class b2ddsj_evt : public TObject{
public:
  b2ddsj_evt(void);
  void ClearEvt(void);
/// exp, run, evtn
  EvtInfo   info;
//  DeltaEMbc dembc;
  double de;
  double mbc;
  double costhBcms;
// SFKW moments and thrust
  EvtShape  shape;
/// B meson decay mode
/// 000 : B0 -> D-       Ds+
/// 001 : B0 -> D-       Dsj+
/// 010 : B0 -> D*-      Ds+
/// 011 : B0 -> D*-      Dsj+
/// 100 : B+ -> anti-D0  Ds+
/// 101 : B+ -> anti-D0  Dsj+
/// 110 : B+ -> anti-D*0 Ds+
/// 111 : B+ -> anti-D*0 Dsj+
  int mode_b;
  bool DsjFl(void) const;
  bool DstFl(void) const;
  bool BchFl(void) const;
/// D meson decay mode
/// 00 : D0 -> K- pi+
/// 01 : D0 -> K- pi+ pi0
/// 02 : D0 -> K- 2pi+ pi-
/// 10 : D+ -> K- pi+ pi+
  int mode_d;
/// Ds meson decay mode
/// 0 : Ds+ -> phi pi+
/// 1 : Ds+ -> anti-K*0 K+, anti-K*0 -> K- pi+
/// 2 : Ds+ -> Ks0 K+
  int mode_ds;
/// D*sj meson decay mode
/// 00 : D*sj+ -> Ds+  gamma
/// 01 : D*sj+ -> Ds+  pi0
/// 02 : D*sj+ -> Ds+  pi+ pi-
/// 10 : D*sj+ -> D*s+ gamma
/// 11 : D*sj+ -> D*s+ pi0
/// 12 : D*sj+ -> D*s+ pi+ pi-
  int mode_dsj;
  bool DsstFl(void) const;
  bool DsjPi0Fl(void) const;
  bool DsjGammaFl(void) const;
/// B meson flavor
///  1 :      B0 or B+
/// -1 : anti-B0 or B-
  int flvb;

/// D0 or D+ mass
  double md;
/// D*sj mass
  double mdsj;
/// D*0 or D*+ mass difference
  double dmdst;
/// Ds+ mass
  double mds;
/// D*s+ mass difference
  double dmdsst;
/// phi or K*0 mass
  double mvec;
/// phi or K*0 momentum components
  double pvec[3];

/// phi or K*0 momentum
  double Pvec(void) const;
/// pt of phi or K*0
  double ptvec(void) const;
/// polar angle of phi or K*0
  double costhvec(void) const;

/// Helicity angle for Dsj -> Ds X
  double cos_hel_dsj;
/// Helicity angle for phi of K*0
  double cos_hel_vec;
/// Helicity angle for V
  double cos_hel_b;

/////////////////////////////
// * D*(s)(j) mesons FSP * //
/////////////////////////////

/// Photon from D*sj -> D(*)s gamma
  GammaInfo gam_dsj;
/// Photon from D*s -> Ds gamma
  GammaInfo gam_dsst;
/// pi0 from D*0 -> D0 pi0
  Pi0Info pi0_dst;
/// pi0 from D*sj -> D(*)s pi0
  Pi0Info pi0_dsj;
/// pi+ from D*- -> D0 pi-
  TrackInfo pi_dst;
/// pi+ from D*sj -> D(*)s pi+ pi-
  TrackInfo pip_dsj;
/// pi- from D*sj -> D(*)s pi+ pi-
  TrackInfo pim_dsj;

/////////////////////
// * D meson FSP * //
/////////////////////

/// K- from D0 -> K- pi+ or D+ -> K- pi+ pi-
  TrackInfo k_d;
/// pi+ from D0 -> K- pi+ or D+ -> K- pi+ pi-
  TrackInfo pi_d;
/// pi- from D+ -> K- 2pi+
  TrackInfo pi2_d;

//////////////////////
// * Ds meson FSP * //
//////////////////////

/// K+ or pi+ from Ds+ -> phi pi+ or Ds+ -> Ks0 K+, K*0 K+
  TrackInfo h_ds;
/// Ks0 from Ds+ -> Ks0 K+
  Ks0Info ks_ds;

/////////////////////////
// * phi and K*0 FSP * //
/////////////////////////

/// h1 from phi -> h1 h2 or K*0 -> h1 h2
  TrackInfo h1_vec;
/// h2 from phi -> h1 h2 or K*0 -> h1 h2
  TrackInfo h2_vec;
/// Beam interaction point and CMS boost vector
  IPBoost ipbst;

  ClassDef(b2ddsj_evt,1)
};

#ifdef __MAKECINT__
#pragma link C++ class b2ddsj_evt;
#endif

#endif

