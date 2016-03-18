#include "b2ddsj_mcevt.h"
#include <cmath>

//#if defined(BELLE_NAMESPACE)
//namespace Belle {
//#endif

ClassImp(b2ddsj_mcevt)

b2ddsj_mcevt::b2ddsj_mcevt(void){
  ClearEvt();
}

void b2ddsj_mcevt::ClearEvt(void){
  genhep.Clear();

  gam_dsj_gen.Clear();
  gam_dsst_gen.Clear();
  pi0_dst_gen.Clear();
  pi0_dsj_gen.Clear();
  pi_dst_gen.Clear();
  pip_dsj_gen.Clear();
  pim_dsj_gen.Clear();

  k_d_gen.Clear();
  pi_d_gen.Clear();
  pi2_d_gen.Clear();

  h_ds_gen.Clear();
  ks_ds_gen.Clear();

  h1_vec_gen.Clear();
  h2_vec_gen.Clear(); 
}

bool b2ddsj_mcevt::Sig(void) const{
  return (b_gen.flag == 1 || b_gen.flag == 5 || b_gen.flag == 10);
}
/// Combinatorial bkg
bool b2ddsj_mcevt::Cmb(void) const{
  return b_gen.flag == -1;
}
/// Non-combinatorial bkg
bool b2ddsj_mcevt::Bkg(void) const{
  return !Sig() && !Cmb();
}

//#if defined(BELLE_NAMESPACE)
//}
//#endif

