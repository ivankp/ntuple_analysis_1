#ifndef ANY_READER_HH
#define ANY_READER_HH

#include <tuple>
#include <algorithm>
#include <cstring>

#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/seq/for_each.hpp>

#include <TLeaf.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include "ivanp/pack.hh"
#include "ivanp/error.hh"
#include "ivanp/debug/type_str.hh"

#define ROOT_LEAF_TYPES \
  (Char_t)(UChar_t)(Short_t)(UShort_t)(Int_t)(UInt_t)(Float_t)(Double_t)(Long64_t)(ULong64_t)(Bool_t)

template <typename> constexpr const char* root_type_str();

#define ROOT_TYPE_STR(r, data, T) \
  template <> \
  constexpr const char* root_type_str<T>() { return BOOST_PP_STRINGIZE(T); }

BOOST_PP_SEQ_FOR_EACH(ROOT_TYPE_STR,,ROOT_LEAF_TYPES)

template <typename... Ts>
class any_reader {
  size_t index;

  template <size_t I>
  using type = ivanp::nth_type<I,Ts...>;

  using common_type = std::common_type_t<Ts...>;

  template <typename T>
  using reader_type = TTreeReaderValue<T>;

  char data[ std::max(sizeof(reader_type<Ts>)...) ];

  template <size_t I=0>
  inline std::enable_if_t<(I<sizeof...(Ts)),size_t>
  get_index(const char* type_name) {
    if (!strcmp(root_type_str<type<I>>(),type_name)) return I;
    else return get_index<I+1>(type_name);
  }
  template <size_t I=0>
  inline std::enable_if_t<(I==sizeof...(Ts)),size_t>
  get_index [[noreturn]] (const char* type_name) {
    throw ivanp::error(ivanp::type_str<any_reader<Ts...>>(),
      " cannot read ",type_name);
  }

  template <size_t I>
  inline auto cast() { return reinterpret_cast<reader_type<type<I>>*>(data); }

  template <typename F, size_t I=0>
  inline std::enable_if_t<(sizeof...(Ts)-I>1)> call(F&& f) {
    if (index==I) f(cast<I>());
    else call<F,I+1>(std::forward<F>(f));
  }
  template <typename F, size_t I=0>
  inline std::enable_if_t<(sizeof...(Ts)-I==1)> call(F&& f) {
    f(cast<I>());
  }

public:
  any_reader(TTreeReader& reader, const char* branch_name)
  : index(get_index(
      reader.GetTree()->GetLeaf(branch_name)->GetTypeName()
    ))
  {
    call([&](auto* p) {
      using T = ivanp::decay_ptr_t<decltype(p)>;
      new(p) T(reader,branch_name);
    });
  }
  ~any_reader() {
    call([](auto* p) {
      using T = ivanp::decay_ptr_t<decltype(p)>;
      p->~T();
    });
  }

  inline common_type operator*() {
    common_type x;
    call([&](auto* p){ x = **p; });
    return x;
  }

  inline const char* GetBranchName() {
    const char* x;
    call([&](auto* p){ x = p->GetBranchName(); });
    return x;
  }
};

using any_float_reader = any_reader<double,float>;

#endif
