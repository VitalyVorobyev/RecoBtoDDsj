#ifndef __B2DDSJ_MCEVT_H__
#define __B2DDSJ_MCEVT_H__

//#include "belle.h"
//#include "mytools/dataStructures.h"
#include "TObject.h"
#include "mytools/datastructures.h"

//#if defined(BELLE_NAMESPACE)
//namespace Belle {
//#endif

class b2ddsj_mcevt : public TObject{
public:
  b2ddsj_mcevt(void);
  void ClearEvt(void);

/// Signal event
  bool Sig(void) const;
/// Combinatorial bkg
  bool Cmb(void) const;
/// Non-combinatorial bkg
  bool Bkg(void) const;

  GenHepEvt genhep;

  GenParticleInfo dsj_gen;
  GenParticleInfo ds_gen;
  GenParticleInfo d_gen;
  GenParticleInfo b_gen;

  GenParticleInfo gam_dsj_gen;
  GenParticleInfo gam_dsst_gen;
  GenParticleInfo pi0_dst_gen;
  GenParticleInfo pi0_dsj_gen;
  GenParticleInfo pi_dst_gen;
  GenParticleInfo pip_dsj_gen;
  GenParticleInfo pim_dsj_gen;

  GenParticleInfo k_d_gen;
  GenParticleInfo pi_d_gen;
  GenParticleInfo pi2_d_gen;

  GenParticleInfo h_ds_gen;
  GenParticleInfo ks_ds_gen;

  GenParticleInfo h1_vec_gen;
  GenParticleInfo h2_vec_gen;

  ClassDef(b2ddsj_mcevt,1)
};

#ifdef __MAKECINT__
#pragma link C++ class b2ddsj_mcevt;
#endif

//#if defined(BELLE_NAMESPACE)
//}
//#endif
#endif

