#include "analyses/diphoton.hh"

class diphoton_analysis: diphoton_analysis_base {
  re_axes ra(bfname);

public:
  diphoton_analysis(analysis_args&& args)
  : diphoton_analysis_base(std::move(args)), ra(conf["binning"]) { }

  void fill_hists();


};
ANALYSIS_FACTORY(diphoton_analysis)

void diphoton_analysis::fill_hists() {

}
