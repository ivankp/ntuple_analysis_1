#ifndef ANALYSIS
#define ANALYSIS "loopsim_analysis.hh"
#include "analysis.hh"

#elif defined(ANALYSIS_GLOBAL) // ===================================

#include "loopsim/LoopSim.hh"

#include ANALYSIS_LS

#elif defined(ANALYSIS_INIT) // =====================================

#include ANALYSIS_LS

const double loopsim_R = conf.at("/jets/loopsim_R"_jp);

branch_reader<Char_t[]> _part(reader,"part");

#elif defined(ANALYSIS_LOOP) // =====================================

// Reset ------------------------------------------------------------
const unsigned np = *_nparticle;
particles.clear();

// Keep track of multi-entry events ---------------------------------
++num_entries;
nlo_bin_t::current_id = *_id;
const bool new_id = (prev_id != nlo_bin_t::current_id);
if (new_id) {
  prev_id = nlo_bin_t::current_id;
  ncount_total += (_ncount ? **_ncount : 1);
  ++num_events;
}

const char part = _part[0];
int nloops = 1, event_order = 1;
switch (part) {
  case 'B': { nloops = 0; event_order = 0; break; }
  case 'R': {
    if (np==njets_min) { nloops = 0; } break;
    if (njets_min > np) throw ivanp::error(
      "njets_min(",njets_min,") > np(",np,')');
    if (njets_min+1 < np) throw ivanp::error(
      "njets_min(",njets_min,")+1 < np(",np,')');
  }
  default: ;
}

loopsim::Event ls_event;
for (unsigned i=0; i<np; ++i) {
  ls_event.particles.emplace_back(_px[i],_py[i],_pz[i],_E[i],[](int kf){
    if (kf==21 or abs(kf)<=6) return 81;
    if (kf==22 or kf==25) return kf;
    throw ivanp::error("unexpected flavour ",kf);
  }(_kf[i]));
}
ls_event.weight = 1;

loopsim::LoopSim _loopsim(
  event_order, nloops, ls_event, loopsim_R, njets_min
  // opt_loopsim_nborn
);

TEST(ent)
unsigned ls_event_i = 0;
while (_loopsim.there_is_a_next_event()) {
  TEST((ls_event_i++))
  const auto& new_ls_event = _loopsim.extract_next_event();

  for (unsigned i=weights.size(); i--; ) // set weights
    nlo_bin_t::weights[i] = weights[i]*new_ls_event.weight;

  // Read particles ---------------------------------------------------
  unsigned n22 = 0; // number of photons
  unsigned n25 = 0; // number of Higgs
  for (const auto& lsp : new_ls_event.particles) {
    const auto kf = lsp.flavour().flavour();
    TEST(kf)
    if (kf == 25) {
      if (n25) throw ivanp::error("more than one Higgs");
      higgs = { lsp[0], lsp[1], lsp[2], lsp[3] };
      ++n25;
    } else if (kf == 22) {
      if (n22>1) throw ivanp::error("more than two photons");
      photons[n22] = { lsp[0], lsp[1], lsp[2], lsp[3] };
      ++n22;
    } else {
      particles.emplace_back(lsp);
    }
  }
  if (!n25 && n22!=2) throw ivanp::error("missing Higgs or photons");
// ------------------------------------------------------------------

#include ANALYSIS_LS
}

#elif defined(ANALYSIS_END) // ======================================

#include ANALYSIS_LS

#endif

