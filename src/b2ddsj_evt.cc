#include "b2ddsj_evt.h"
#include <cmath>

//#if defined(BELLE_NAMESPACE)
//namespace Belle {
//#endif

ClassImp(b2ddsj_evt)

b2ddsj_evt::b2ddsj_evt(void){
  ClearEvt();
}

void b2ddsj_evt::ClearEvt(void){
  info.Clear();
  de  = -2.;
  mbc = 0.;
  costhBcms = -2.;
  shape.Clear();
  mode_b = -1; mode_d = -1; mode_ds = -1; mode_dsj = -1; flvb = 0;
  md = 0; mdsj = 0; dmdst = 0; mds = 0; dmdsst = 0;
  mvec = 0;
  pvec[0] = -99; pvec[1] = -99; pvec[2] = -99;
  cos_hel_dsj = -2.; cos_hel_vec = -2.; cos_hel_b = -2.;

  gam_dsj.Clear();
  gam_dsst.Clear();
  pi0_dst.Clear();
  pi0_dsj.Clear();
  pi_dst.Clear();
  pip_dsj.Clear();
  pim_dsj.Clear();

  k_d.Clear();
  pi_d.Clear();
  pi2_d.Clear();

  h_ds.Clear();
  ks_ds.Clear();

  h1_vec.Clear();
  h2_vec.Clear();

  ipbst.Clear();
}

bool b2ddsj_evt::DsjFl(void) const{
  if(mode_b < 0) return false;
  return mode_b % 2;
}

bool b2ddsj_evt::DstFl(void) const{
  if(mode_b < 0) return false;
  return (mode_b % 100) / 10;
}

bool b2ddsj_evt::BchFl(void) const{
  if(mode_b < 0) return false;
  return mode_b / 100;
}

bool b2ddsj_evt::DsstFl(void) const{
  if(mode_dsj < 0) return false;
  return mode_dsj / 10;
}

bool b2ddsj_evt::DsjPi0Fl(void) const{
  if(mode_dsj < 0) return false;
  return (mode_dsj % 10) == 1;
}

bool b2ddsj_evt::DsjGammaFl(void) const{
  if(mode_dsj < 0) return false;
  return (mode_dsj % 10) == 0;
}

double b2ddsj_evt::Pvec(void) const{
  return sqrt(pvec[0]*pvec[0]+pvec[1]*pvec[1]+pvec[2]*pvec[2]);
}

double b2ddsj_evt::ptvec(void) const{
  return sqrt(pvec[0]*pvec[0]+pvec[1]*pvec[1]);
}

double b2ddsj_evt::costhvec(void) const{
  return pvec[3]/Pvec();
}

//#if defined(BELLE_NAMESPACE)
//}
//#endif

