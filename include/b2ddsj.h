#include "particle/Particle.h"
#include "particle/utility.h"
#include "event/BelleEvent.h"
#include "basf/module.h"
#include "basf/module_descr.h"
#include "ip/IpProfile.h"
#include "particle/combination.h"

#include "benergy/BeamEnergy.h"

//#include "hamlet/Hamlet.h"

#include <vector>
#include "belle.h"

#include "TTree.h"
#include "TFile.h"

#include "b2ddsj_evt.h"
#include "b2ddsj_mcevt.h"
#include "b2ddsj_evt_dict.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class b2ddsj : public Module{
public:
  b2ddsj(void):
  clock(0),ntuple_flag(1),n_good_b(0) {}
  ~b2ddsj(void) {}
  void init(int *);
  void term(void);
  void disp_stat(const char*){}
  void hist_def(void) {}
  void begin_run(BelleEvent*, int*);
  void event(BelleEvent*, int*);
  void end_run(BelleEvent*, int*){}

  int m_mode;//0 -> Data
             //1 -> Signal MC
             //2 -> Genegic M
  int ntuple_flag;
  char ofile[1024];

private:
  int PrepareDs(void);
  int PrepareD0(void);
  int PrepareDp(void);
  int PrepareDst(void);
  int PrepareDsj(void);
  int PrepareB(void);

  b2ddsj_evt*   tevt;
  b2ddsj_mcevt* tmcevt;
  void SetupTree(void);

  bool IsMC(void);
  bool mc_flag;
  TTree* TEvent;
  TFile* tfile;
  int clock;
  int n_good_b;
  void begin_event(BelleEvent* evptr);

  std::vector<Particle> pipl,piml,kpl,kml;
  std::vector<Particle> ks0l;
  std::vector<Particle> phil;
  std::vector<Particle> pi0l;
  std::vector<Particle> gammal;
  std::vector<Particle> d0l,d0bl;
  std::vector<Particle> dpl,dml;
  std::vector<Particle> dst0l,dst0bl;
  std::vector<Particle> dstpl,dstml;
  std::vector<Particle> kst0l,kst0bl;
  std::vector<Particle> dspl,dsml;
  std::vector<Particle> dsstpl,dsstml;
  std::vector<Particle> dsjpl,dsjml;
  std::vector<Particle> bl;

  int DumpEvent(void);
  void FillEvt(Particle& b);
  void FillVec(const Particle& vec);
};

#if defined(BELLE_NAMESPACE)
}
#endif

