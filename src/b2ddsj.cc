#include "b2ddsj.h"

#include "mytools/recotools.h"
#include "mytools/uisetter.h"
#include "mytools/combinator.h"
#include "mytools/geninfo.h"

#include <iostream>

using namespace std;

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

const bool dump = false;
extern "C" Module_descr *mdcl_b2ddsj()
{
  b2ddsj *module = new b2ddsj;
  Module_descr *dscr = new Module_descr("b2ddsj",module);
  BeamEnergy::define_global(dscr);
  dscr->define_param("mode","mode",&module->m_mode);
  dscr->define_param("ofile","ofile",1024,module->ofile);
  dscr->define_param("ntuple_flag","ntuple_flag",&module->ntuple_flag);
  return dscr;
}

void b2ddsj::init(int *){
  std::cout << "init" << endl;
  Ptype dummy("VPHO");
//  if(ntuple_flag) Hamlet::init();
  return;
}

bool b2ddsj::IsMC(void){
  Belle_runhead_Manager &rhd_mgr = Belle_runhead_Manager::get_manager();
  std::vector<Belle_runhead>::const_iterator rhd = rhd_mgr.begin();
  if(rhd == rhd_mgr.end()){
    fprintf(stderr,"Constructor: Cannot access to Belle_runhead\n");
    return false;
  } else{
    if(rhd->ExpMC() == 1) return false;//real data
    else                  return true;
  }
  return false;
}

void b2ddsj::begin_run(BelleEvent* evptr, int *status){
  IpProfile::begin_run();
  BeamEnergy::begin_run();
  eid::init_data();
//  if(ntuple_flag){
///    Hamlet::begin_run(Hamlet::MULT_DIM_LH);
//  } else{
//    UISetter::m_rphi_svd = 0;// no SVD hits requirement at skim level
//    UISetter::m_rz_svd   = 0;
//  }
  mc_flag = IsMC();
  if(ntuple_flag && !clock){
    clock = 1;
    tfile = new TFile(ofile,"RECREATE");
    cout << ofile << endl;
    SetupTree();
  }
  Combinator::SetNT(ntuple_flag);
  Combinator::SetMC(mc_flag);
  UISetter::SetMC(mc_flag);

  return;
}

void b2ddsj::SetupTree(void){
  TEvent = new TTree("tevt","tevt");
  tevt = new b2ddsj_evt();
  TEvent->Branch("evt",tevt,"evt/b2ddsj_evt");
  if(m_mode){
    tmcevt = new b2ddsj_mcevt();
    TEvent->Branch("mcevt",tmcevt,"mcevt/b2ddsj_mcevt");
  }
  cout << "TTree tevt is initialized" << endl;
}

void b2ddsj::term(void){
  cout << "Term " << ofile << endl;
  if(ntuple_flag){
    cout << "Term TEvent: " << TEvent->GetEntries() << endl;
    TEvent->Write();
    tfile->Close();
  }
  cout << "# of Good B0: " << n_good_b << endl;
}

void b2ddsj::begin_event(BelleEvent* evptr){
  Belle_event_Manager& evtmgr = Belle_event_Manager::get_manager();
  Belle_event& evthead = *evtmgr.begin();
  const int evtn = evthead.EvtNo() & 0xfffffff;
  if(!(evtn%100)) cout << "run: " << evthead.RunNo() << ", evtn: " << evtn << endl;
  if(!IpProfile::usable()){ std::cout<<"IP profile is not available!!!\n"; return; }
}

void b2ddsj::event(BelleEvent* evptr, int *status){
  *status = 0; begin_event(evptr);
  tevt->ClearEvt();
  if(m_mode) tmcevt->ClearEvt();
  const bool dump = false;
  Combinator::make_kpi(pipl,piml,kpl,kml);
  Combinator::make_ks(ks0l);
  Combinator::make_phitokk(phil,kpl,kml);
  Combinator::make_kstar0(kst0l,kst0bl,kpl,piml,kml,pipl);


  if(!PrepareDs()) return;
  const int nD0  = PrepareD0();
  const int nDpm = PrepareDp();
  if(!(nD0+nDpm)) return;
  if(!nD0 && !(dspl.size()+dml.size()) && !(dsml.size()+dpl.size())) return;

  Combinator::make_pi0(pi0l);
  Combinator::make_gamma(gammal);

  PrepareDst();
  if(!PrepareDsj()) return;

  const int nB = PrepareB();
  if(!nB) return;
  n_good_b += nB;
  if(dump) DumpEvent();
// ** ** //
  if(ntuple_flag){
    for(unsigned i=0; i<bl.size(); i++){
      FillEvt(bl[i]);
      TEvent->Fill();
    }
  } else{
    *status = 1;
  }

  return;
}

int b2ddsj::DumpEvent(void){
  cout << "Ds+: " << dspl.size() << ", Ds-: " << dsml.size() << endl;
  cout << "D0 : " << d0l.size() + d0bl.size() << ", D+ : " << dpl.size() << ",D-: " << dml.size() << endl;
  cout << "pi0: " << pi0l.size() << ", gamma: " << gammal.size() << endl;
  cout << "     D*0: " << dst0l.size()  << endl;
  cout << "anti-D*0: " << dst0bl.size() << endl;
  cout << "     D*+: " << dstpl.size()  << endl;
  cout << "     D*-: " << dstml.size()  << endl;
  cout << "    Ds*+: " << dsstpl.size() << endl;
  cout << "    Ds*-: " << dsstml.size() << endl;
  cout << "   Dsj*+: " << dsjpl.size() << endl;
  cout << "   Dsj*-: " << dsjml.size() << endl;
  cout << "B: " << bl.size() << endl;
  return 0;
}

int b2ddsj::PrepareDst(void){
  Combinator::make_dsstar(dsstpl,dspl,gammal);
  Combinator::make_dsstar(dsstml,dsml,gammal);

  Combinator::make_dstar(dstpl,d0l, pipl);
  Combinator::make_dstar(dstml,d0bl,piml);

  Combinator::make_dstar(dst0l, d0l, pi0l);
  Combinator::make_dstar(dst0bl,d0bl,pi0l);
  return 0;
}

int b2ddsj::PrepareDs(void){
  // Ds+ -> phi pi+ (mode 0)
  Combinator::make_dstoh0hp(dspl,dsml,phil,pipl,piml,0);
  // Ds+ -> K*0 K+  (mode 1)
  vector<Particle> dspl_Kst0K,dsml_Kst0K;
  Combinator::make_dstoKK(dspl_Kst0K,dsml_Kst0K,kst0bl,kpl,kst0l,kml,1);
  dspl.insert(dspl.end(),dspl_Kst0K.begin(),dspl_Kst0K.end());
  dsml.insert(dsml.end(),dsml_Kst0K.begin(),dsml_Kst0K.end());
  // Ds+ -> Ks0 K+  (mode 2)
  vector<Particle> dspl_Ks0K,dsml_Ks0K;
  Combinator::make_dstoh0hp(dspl_Ks0K,dsml_Ks0K,ks0l,kpl,kml,2);
  dspl.insert(dspl.end(),dspl_Ks0K.begin(),dspl_Ks0K.end());
  dsml.insert(dsml.end(),dsml_Ks0K.begin(),dsml_Ks0K.end());
  return dspl.size() + dsml.size();
}

int b2ddsj::PrepareD0(void){
  // D0 -> K- pi+
  Combinator::make_d0tokpi(d0l,d0bl,kml,pipl,kpl,piml);
  return d0l.size() + d0bl.size();
}

int b2ddsj::PrepareDp(void){
  // Dp -> K- pi+ pi-
  Combinator::make_dptokpipi(dpl,dml,kml,pipl,kpl,piml);
  return dpl.size() + dml.size();
}

int b2ddsj::PrepareDsj(void){
  const bool dump = false;
  dsjpl.clear(); dsjml.clear();
// Dsj -> Ds gamma
  vector<Particle> dsp_gam_l,dsm_gam_l;
  if(dump) cout << "PrepareDsj1" << endl;
  Combinator::make_dsjtodsx(dsp_gam_l,dsm_gam_l,dspl,dsml,gammal,0);
  dsjpl.insert(dsjpl.end(),dsp_gam_l.begin(),dsp_gam_l.end());
  dsjml.insert(dsjml.end(),dsm_gam_l.begin(),dsm_gam_l.end());
// Dsj -> Ds pi0
  vector<Particle> dsp_pi0_l,dsm_pi0_l;
  if(dump) cout << "PrepareDsj2" << endl;
  Combinator::make_dsjtodsx(dsp_pi0_l,dsm_pi0_l,dspl,dsml,pi0l,1);
  dsjpl.insert(dsjpl.end(),dsp_pi0_l.begin(),dsp_pi0_l.end());
  dsjml.insert(dsjml.end(),dsm_pi0_l.begin(),dsm_pi0_l.end());
// Dsj -> Ds pi+ pi-
  vector<Particle> dsp_pipi_l,dsm_pipi_l;
  if(dump) cout << "PrepareDsj3" << endl;
  Combinator::make_dsjtodsxy(dsp_pi0_l,dsm_pi0_l,dspl,dsml,pipl,piml,2);
  dsjpl.insert(dsjpl.end(),dsp_pipi_l.begin(),dsp_pipi_l.end());
  dsjml.insert(dsjml.end(),dsm_pipi_l.begin(),dsm_pipi_l.end());
// Dsj -> Ds* gamma
  vector<Particle> dsstp_gam_l,dsstm_gam_l;
  if(dump) cout << "PrepareDsj4" << endl;
  Combinator::make_dsjtodsx(dsstp_gam_l,dsstm_gam_l,dsstpl,dsstml,gammal,10);
  dsjpl.insert(dsjpl.end(),dsstp_gam_l.begin(),dsstp_gam_l.end());
  dsjml.insert(dsjml.end(),dsstm_gam_l.begin(),dsstm_gam_l.end());
// Dsj -> Ds* pi0
  vector<Particle> dsstp_pi0_l, dsstm_pi0_l;
  if(dump) cout << "PrepareDsj5" << endl;
  Combinator::make_dsjtodsx(dsstp_pi0_l,dsstm_pi0_l,dsstpl,dsstml,pi0l,11);
  dsjpl.insert(dsjpl.end(),dsstp_pi0_l.begin(),dsstp_pi0_l.end());
  dsjml.insert(dsjml.end(),dsstm_pi0_l.begin(),dsstm_pi0_l.end());
// Dsj -> Ds* pi+ pi-
  vector<Particle> dsstp_pipi_l, dsstm_pipi_l;
  if(dump) cout << "PrepareDsj6" << endl;
  Combinator::make_dsjtodsxy(dsstp_pi0_l,dsstm_pi0_l,dsstpl,dsstml,pipl,piml,12);
  dsjpl.insert(dsjpl.end(),dsstp_pipi_l.begin(),dsstp_pipi_l.end());
  dsjml.insert(dsjml.end(),dsstm_pipi_l.begin(),dsstm_pipi_l.end());
  return dsjpl.size() + dsjml.size();
}

int b2ddsj::PrepareB(void){
  const bool dump = false;
  bl.clear();
  // B0 -> D- Ds+ and c.c.
  vector<Particle> b0ddsl;
  if(dump) cout << "PrepareB1" << endl;
  Combinator::make_b0toxy(b0ddsl,dml,dspl,dpl,dsml,0);
  bl.insert(bl.end(),b0ddsl.begin(),b0ddsl.end());
  // B0 -> D- D*sj+ and c.c.
  vector<Particle> b0ddsjl;
  if(dump) cout << "PrepareB2" << endl;
  Combinator::make_b0toxy(b0ddsjl,dml,dsjpl,dpl,dsjml,1);
  bl.insert(bl.end(),b0ddsjl.begin(),b0ddsjl.end());
  // B0 -> D*- D*sj+ and c.c
  vector<Particle> b0dstdsjl;
  if(dump) cout << "PrepareB3" << endl;
  Combinator::make_b0toxy(b0ddsjl,dstml,dsjpl,dstpl,dsjml,11);
  bl.insert(bl.end(),b0dstdsjl.begin(),b0dstdsjl.end());
  // B+ -> anti-D0 Ds+ and c.c
  vector<Particle> bpd0dsl,bmd0dsl;
  if(dump) cout << "PrepareB4" << endl;
  cout << d0bl.size() << " " << dspl.size() << " " << d0l.size() << " " << dsml.size() << endl;
  Combinator::make_bptoxy(bpd0dsl,d0bl,dspl,100);
  Combinator::make_bptoxy(bmd0dsl,d0l, dsml,100);
  bl.insert(bl.end(),bpd0dsl.begin(),bpd0dsl.end());
  bl.insert(bl.end(),bmd0dsl.begin(),bmd0dsl.end());
  // B+ -> anti-D0 D*sj+ and c.c
  vector<Particle> bpd0dsjl,bmd0dsjl;
  if(dump) cout << "PrepareB5" << endl;
  Combinator::make_bptoxy(bpd0dsjl,d0bl,dsjpl,101);
  Combinator::make_bptoxy(bmd0dsjl,d0l, dsjml,101);
  bl.insert(bl.end(),bpd0dsjl.begin(),bpd0dsjl.end());
  bl.insert(bl.end(),bmd0dsjl.begin(),bmd0dsjl.end());
  // B+ -> anti-D*0 D*sj+ and c.c
  vector<Particle> bpdst0dsjl,bmdst0dsjl;
  if(dump) cout << "PrepareB6" << endl;
  Combinator::make_bptoxy(bpd0dsjl,dst0bl,dsjpl,111);
  Combinator::make_bptoxy(bmd0dsjl,dst0l, dsjml,111);
  bl.insert(bl.end(),bpd0dsjl.begin(),bpd0dsjl.end());
  bl.insert(bl.end(),bmd0dsjl.begin(),bmd0dsjl.end());
  return bl.size();
}

void b2ddsj::FillEvt(Particle& b){
  RTools::FillEvtInfo(tevt->info);
  RTools::FillIPBoost(tevt->ipbst);
  RTools::FillShape(b,tevt->shape);

  if(m_mode){// MC
    RTools::FillGenHepEvt(tmcevt->genhep);
    setMCtruth(b);
    RTools::FillGenPInfo(b,tmcevt->b_gen);
  }

  const B0UserInfo& binfo = (B0UserInfo&)b.userInfo();
  tevt->mbc       = binfo.Mbc();
  tevt->de        = binfo.deltaE();
  tevt->mode_b    = binfo.Mode();
  tevt->costhBcms = RTools::cosThetaCMS(b);

  if(tevt->DsjFl()){
    const Particle&  Dsj = b.child(1);
    cout << "Dsj: " << Dsj.pType().lund() << endl;
    if(m_mode) RTools::FillGenPInfo(Dsj,tmcevt->dsj_gen);
    const DUserInfo& dsjinfo = (DUserInfo&)(Dsj.userInfo());
    tevt->mode_dsj    = dsjinfo.Mode();
    tevt->mdsj        = Dsj.p().m();
    tevt->cos_hel_dsj = RTools::Helicity(Dsj);
    switch(tevt->mode_dsj % 10){
    case 0:{// Dsj* -> X gamma
      const Particle& gamma = Dsj.child(1);
      RTools::FillGamma(gamma,tevt->gam_dsj);
      if(m_mode) RTools::FillGenPInfo(gamma,tmcevt->gam_dsj_gen);
      break;}
    case 1:{// Dsj* -> X pi0
      const Particle& pi0 = Dsj.child(1);
      RTools::FillPi0(pi0,tevt->pi0_dsj);
      if(m_mode) RTools::FillGenPInfo(pi0,tmcevt->pi0_dsj_gen);
      break;}
    case 2:{// Dsj* -> X pi+ pi-
      const Particle& pip = Dsj.child(1);
      const Particle& pim = Dsj.child(2);
      RTools::FillTrk(pip,tevt->pip_dsj);
      RTools::FillTrk(pim,tevt->pim_dsj);
      if(m_mode){
        RTools::FillGenPInfo(pip,tmcevt->pip_dsj_gen);
        RTools::FillGenPInfo(pim,tmcevt->pim_dsj_gen);
      }
      break;}
    default:
      cout << "b2ddsj::FillEvt: unknown Dsj* mode " << tevt->mode_dsj << endl;
      break;
    }
  }

  if(tevt->DstFl()){
    tevt->dmdst  = b.child(0).p().m() - b.child(0).child(0).p().m();
    if(tevt->BchFl()){
      const Particle& pi0 = b.child(0).child(1);
      RTools::FillPi0(pi0,tevt->pi0_dst);
//      if(m_mode) RTools::FillGenPInfo(pi0,tmcevt->pi0_dst_gen);
    } else{
      const Particle& pi = b.child(0).child(1);
      RTools::FillTrk(pi,tevt->pi_dst);
//      if(m_mode) RTools::FillGenPInfo(pi,tmcevt->pi_dst_gen);
    }
  }
  if(tevt->DsstFl()){
    tevt->dmdsst = b.child(1).p().m() - b.child(1).child(0).p().m();
    const Particle& gamma = b.child(1).child(1);
    RTools::FillGamma(gamma,tevt->gam_dsst);
//    if(m_mode) RTools::FillGenPInfo(gamma,tmcevt->gam_dsst_gen);
  }

  const Particle& D  = tevt->DstFl() ? b.child(0).child(0) : b.child(0);
  const Particle& Ds = tevt->DsjFl() ? (tevt->DsstFl() ? b.child(1).child(0).child(0) : b.child(1).child(0)) : b.child(1);
  const DUserInfo& dinfo  = (DUserInfo&)D.userInfo();
  const DUserInfo& dsinfo = (DUserInfo&)Ds.userInfo();
  tevt->mode_d  = dinfo.Mode();
  tevt->mode_ds = dsinfo.Mode();
  tevt->md  = D.p().m();
  tevt->mds = Ds.p().m();

  switch(tevt->mode_d){
  case 10:// D0 -> K- pi+
    RTools::FillTrk(D.child(0),tevt->k_d);
    RTools::FillTrk(D.child(1),tevt->pi_d);
    break;
  case 110:// D+ -> K- 2pi+
    RTools::FillTrk(D.child(0),tevt->k_d);
    RTools::FillTrk(D.child(1),tevt->pi_d);
    RTools::FillTrk(D.child(2),tevt->pi2_d);
    break;
  default:
    cout << "b2ddsj::FillEvt: unknown D meson decay: " << tevt->mode_d << endl;
    break;
  }

  RTools::FillTrk(Ds.child(1),tevt->h_ds);
  switch(tevt->mode_ds){
  case 0:// Ds+ -> phi pi+
    cout << "Ds+ mode " << tevt->mode_ds << endl;
    FillVec(Ds.child(0));
    break;
  case 1:// Ds+ -> K*0 K+
    FillVec(Ds.child(0));
    break;
  case 2:// Ds+ -> Ks0 K+
    RTools::FillKs0(Ds.child(0),tevt->ks_ds);
    break;
  default:
    cout << "b2ddsj::FillEvt: unknown Ds meson decay: " << tevt->mode_ds << endl;
    break;
  }
}

void b2ddsj::FillVec(const Particle& vec){
  cout << "FillVec: " << vec.pType().lund() << " " << vec.nChildren() << endl;
  RTools::FillTrk(vec.child(0),tevt->h1_vec);
  RTools::FillTrk(vec.child(1),tevt->h2_vec);
  tevt->mvec    = vec.p().m();
  tevt->pvec[0] = vec.p().x();
  tevt->pvec[1] = vec.p().y();
  tevt->pvec[2] = vec.p().z();
  tevt->cos_hel_vec = RTools::Helicity(vec);
}

#if defined(BELLE_NAMESPACE)
}
#endif

