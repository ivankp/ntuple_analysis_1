#include "reweighter.hh"

#include <cmath>
#include <array>
#include <map>
#include <functional>

#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/seq/enum.hpp>
#include <boost/preprocessor/seq/transform.hpp>
#include <boost/preprocessor/seq/variadic_seq_to_seq.hpp>

#include <LHAPDF/LHAPDF.h>

#include "branch_reader.hh"

#define REWEIGHTING_BRANCHES \
  ((Int_t)(nparticle)) \
  ((Double_t[],Float_t)(px)) \
  ((Double_t[],Float_t)(py)) \
  ((Double_t[],Float_t)(pz)) \
  ((Double_t[],Float_t)(E)) \
  ((Int_t[])(kf)) \
  ((Double_t)(alphas)) \
  ((Double_t)(weight2)) \
  ((Double_t)(me_wgt)) \
  ((Double_t)(me_wgt2)) \
  ((Double_t)(x1)) \
  ((Double_t)(x2)) \
  ((Double_t)(x1p)) \
  ((Double_t)(x2p)) \
  ((Int_t)(id1)) \
  ((Int_t)(id2)) \
  ((Double_t)(fac_scale)) \
  ((Double_t)(ren_scale)) \
  ((Double_t[])(usr_wgts)) \
  ((Char_t)(alphasPower)) \
  ((Char_t[])(part))

#define PP_EXPAND(...) __VA_ARGS__
#define STRIP(X) PP_EXPAND( PP_EXPAND X )

#define SEQ_ELEM(I,SEQ) \
  STRIP(BOOST_PP_SEQ_ELEM(I,BOOST_PP_VARIADIC_SEQ_TO_SEQ(SEQ)))

struct branches {
#define MAKE_BRANCH_READER(r, data, x) \
  branch_reader<SEQ_ELEM(0,x)> SEQ_ELEM(1,x);

  BOOST_PP_SEQ_FOR_EACH(MAKE_BRANCH_READER,,REWEIGHTING_BRANCHES)
#undef MAKE_BRANCH_READER

#define MAKE_BRANCH_READER(r, data, x) \
  SEQ_ELEM(1,x)( reader, BOOST_PP_STRINGIZE(SEQ_ELEM(1,x)) )

  branches(TTreeReader& reader)
  : BOOST_PP_SEQ_ENUM(
      BOOST_PP_SEQ_TRANSFORM(MAKE_BRANCH_READER,,REWEIGHTING_BRANCHES))
  { }
#undef MAKE_BRANCH_READER
};

std::vector<std::unique_ptr<LHAPDF::PDF>>
make_pdfs(const std::string& name, bool variations) {
  if (variations) {
    auto pdfs = LHAPDF::mkPDFs(name);
    return { pdfs.begin(), pdfs.end() };
  } else {
    std::vector<std::unique_ptr<LHAPDF::PDF>> pdfs;
    pdfs.emplace_back( LHAPDF::mkPDF(name) );
    return std::move(pdfs);
  }
}

std::map<
  std::string,
  std::function<double(double k, branches& b)>
> scale_functions {
};

struct reweighter_impl: branches {

  reweighter::args_struct args;
  const std::vector<std::unique_ptr<LHAPDF::PDF>> pdfs;
  LHAPDF::PDF *pdf = nullptr;

  decltype(scale_functions)::mapped_type& scale_f;

  std::vector<double> scale_values;
  std::vector<unsigned> Kri, Kfi;

  struct ren_vars_struct { double ar, w0; };
  struct fac_vars_struct { double  m, ff; };
  std::vector<ren_vars_struct> ren_vars;
  std::vector<fac_vars_struct> fac_vars;

  std::vector<double> weights;

  struct fac_calc_struct {
    LHAPDF::PDF *pdf;
    std::array<int,2> id;
    double muF;
    std::array<double,2> x, xp;

    static constexpr std::array<int,10> quarks {
      1,-1, 2,-2, 3,-3, 4,-4, 5,-5
    };

    double fxQ(int id, double x1, double x2) const {
      return pdf->xfxQ(id, x1, muF)/x2;
    }
    double fxQ(int id, double x) const {
      return pdf->xfxQ(id, x, muF)/x;
    }
    double fr(unsigned r) const {
      return fxQ(id[r],x[r]);
    }
    double fr1(unsigned r) const {
      if (id[r] != 21) return fxQ(id[r],x[r]);
      else {
        double f = 0.;
        for (int q : quarks) f += fxQ(q,x[r]);
        return f;
      }
    };
    double fr2(unsigned r) const {
      if (id[r] != 21) return fxQ(id[r],x[r]/xp[r],x[r]);
      else {
        double f = 0.;
        for (int q : quarks) f += fxQ(q,x[r]/xp[r],x[r]);
        return f;
      }
    };
    double fr3(unsigned r) const {
      return fxQ(21,x[r]);
    }
    double fr4(unsigned r) const {
      return fxQ(21,x[r]/xp[r],x[r]);
    }
  } fc;

  void fac_calc(unsigned i) {
    auto& v = fac_vars[i];
    fc.muF = scale_values[i];
    fc.id  = { *id1, *id2 };
    fc.x   = { * x1, * x2 };

    const double f[2] = { fc.fr(0), fc.fr(1) };
    v.ff = f[0]*f[1];

    if (part[0]=='I') {
      fc.xp = { *x1p, *x2p };
      const double lf = 2.*std::log(fc.muF/(*fac_scale));
      double w[8];
      for (int i=0; i<8; ++i) w[i] = usr_wgts[i+2] + usr_wgts[i+10]*lf;

      v.m = ( fc.fr1(0)*w[0] + fc.fr2(0)*w[1]
            + fc.fr3(0)*w[2] + fc.fr4(0)*w[3] )*f[1]
          + ( fc.fr1(1)*w[4] + fc.fr2(1)*w[5]
            + fc.fr3(1)*w[6] + fc.fr4(1)*w[7] )*f[0];
    } else v.m = 0.;
  }

  void ren_calc(unsigned i) {
    auto& v = ren_vars[i];
    const double muR = scale_values[i];

    v.ar = std::pow(pdf->alphasQ(muR)/(*alphas), (*alphasPower));

    if (part[0]=='V' || part[0]=='I') {
      const double lr = 2.*std::log(muR/(*ren_scale));
      v.w0 = (*me_wgt) + lr*usr_wgts[0] + 0.5*lr*lr*usr_wgts[1];
    } else {
      v.w0 = *me_wgt2;
    }
  }

  double combine(unsigned i) { // i is scale index
    const auto& scale = args.Ki[i];
    double w;

    if (scale.fac) {
      const auto& fac = fac_vars[*scale.fac];
      w = fac.m;
      if (scale.ren) {
        const auto& ren = ren_vars[*scale.ren];
        w += ren.w0 * fac.ff;
      } else {
        w += (*me_wgt2) * fac.ff;
      }
    } else {
      w = *weight2;
    }
    if (scale.ren) {
      const auto& ren = ren_vars[*scale.ren];
      w *= ren.ar;
    }

    return w;
  }

  void operator()() {
    for (unsigned i=0; i<args.K.size(); ++i)
      scale_values[i] = scale_f(args.K[i],*this);

    pdf = fc.pdf = pdfs[0].get();
    for (unsigned i : Kri) ren_calc(i);
    for (unsigned i : Kfi) fac_calc(i);

    // scale variations
    for (unsigned i=0; i<args.Ki.size(); ++i)
      weights[i] = combine(i);

    // pdf variations
    for (unsigned i=1; i<pdfs.size(); ++i) {
      pdf = fc.pdf = pdfs[i].get();
      fac_calc(0);
      ren_calc(0);
      weights[args.Ki.size()+i] = combine(0);
    }
  }

  reweighter_impl(TTreeReader& reader, reweighter::args_struct args)
  : branches(reader),
    args(std::move(args)),
    pdfs(make_pdfs(args.pdf,args.pdf_var)),
    scale_f(scale_functions[args.scale]),
    scale_values(args.K.size()), Kri(), Kfi(),
    ren_vars(), fac_vars(),
    weights(args.Ki.size()+pdfs.size()-1)
  {
    Kri.reserve(args.K.size());
    Kfi.reserve(args.K.size());
    for (const auto& ki : args.Ki) {
      if (ki.ren)
        if (std::find(Kri.begin(),Kri.end(),*ki.ren)!=Kri.end())
          Kri.push_back(*ki.ren);
      if (ki.fac)
        if (std::find(Kfi.begin(),Kfi.end(),*ki.ren)!=Kfi.end())
          Kfi.push_back(*ki.fac);
    }
    ren_vars.resize(Kri.size());
    fac_vars.resize(Kfi.size());
  }
};

reweighter::reweighter(TTreeReader& reader, args_struct args)
: impl(new reweighter_impl(reader, std::move(args))) { }
reweighter::~reweighter() { delete impl; }

void reweighter::args_struct::add_scale(const ren_fac<double>& k) {
  ren_fac<unsigned> ki { };
  if (k.ren) {
    auto it = std::find(K.begin(), K.end(), *k.ren);
    if (it==K.end()) {
      ki.ren = K.size();
      K.push_back(*k.ren);
    } else {
      ki.ren = it - K.begin();
    }
  }
  if (k.fac) {
    auto it = std::find(K.begin(), K.end(), *k.fac);
    if (it==K.end()) {
      ki.fac = K.size();
      K.push_back(*k.fac);
    } else {
      ki.fac = it - K.begin();
    }
  }
  Ki.push_back(ki);
}

unsigned reweighter::nweights() const {
  return impl->weights.size();
}
double reweighter::operator[](unsigned i) const {
  return impl->weights[i];
}
std::string reweighter::name(unsigned i) const {
  const auto *pdf = impl->pdfs[i].get();
  return pdf->set().name() + "_" + std::to_string(pdf->memberID());
}

