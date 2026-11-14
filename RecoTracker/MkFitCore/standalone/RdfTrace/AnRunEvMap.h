#ifndef RecoTracker_MkFitCore_standalone_RdfTrace_AnRunEvMap_h
#define RecoTracker_MkFitCore_standalone_RdfTrace_AnRunEvMap_h

// Projections of an Event's member vectors into RDataFrame columns.
//
// One RDF entry is one mkfit::Event, so every per-object quantity reaches the
// dataframe as an RVec built from one of the Event's vectors -- trSeeds_,
// trCandMetas_, simTracks_, any of them. These are the helpers and the macros
// that do it. They lived in AnRun.cc until a second translation unit needed
// them; nothing about them is specific to one analysis.

#include "RecoTracker/MkFitCore/standalone/Event.h"
#include "ROOT/RVec.hxx"

#include <cassert>
#include <vector>

//==============================================================================
#pragma region Map/Gather/...
//==============================================================================

namespace { // map, gather, compress

  // using RVecI = ROOT::RVec<int>;

  // ---------------------------------------------------------------------------
  // map_with_member:
  // ---------------------------------------------------------------------------
  template <typename VEC_T, typename F>
  auto map_with_member(const VEC_T& v, F f) {
    using T = std::decay_t<decltype(f(*v.begin()))>;
    ROOT::RVec<T> out;
    out.reserve(v.size());
    for (auto& x : v) out.push_back(f(x));
    return out;
  }
  // ---------------------------------------------------------------------------
  // map_with_func:
  // ---------------------------------------------------------------------------
  template <typename T, typename F>
  auto map_with_func(const std::vector<T>& v, const mkfit::Event* ev, F func) {
    using R = std::decay_t<decltype(func(ev, std::declval<T>()))>;
    ROOT::RVec<R> out;
    out.reserve(v.size());
    for (const auto& elem : v) out.push_back(func(ev, elem));
    return out;
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  template <typename VEC_T, typename F>
  auto gather_with_member(const VEC_T& v, F f, const ROOT::RVec<int>& idx_vec) {
    using T = std::decay_t<decltype(f(*v.begin()))>;
    ROOT::RVec<T> out;
    out.reserve(idx_vec.size());
    for (auto &i : idx_vec) out.push_back(f(v[i]));
    return out;
  }
  // ---------------------------------------------------------------------------
  // gather_with_func: Call function on selected indices only
  // ---------------------------------------------------------------------------
  template <typename T, typename F>
  auto gather_with_func(const std::vector<T>& v, const mkfit::Event* ev,
                        F func, const ROOT::RVec<int>& idx_vec) {
    using R = std::decay_t<decltype(func(ev, std::declval<T>()))>;
    ROOT::RVec<R> out;
    out.reserve(idx_vec.size());
    for (int i : idx_vec) out.push_back(func(ev, v[i]));
    return out;
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  template <typename VEC_T, typename F>
  auto compress_with_member(const VEC_T& v, F f, const ROOT::RVec<int>& mask) {
    assert(v.size() == mask.size());
    using T = std::decay_t<decltype(f(*v.begin()))>;
    ROOT::RVec<T> out;
    const size_t n = v.size();
    out.reserve(n / 16);
    for (size_t i = 0; i < n; ++i) {
      if (mask[i]) out.push_back(f(v[i]));
    }
    return out;
  }
  // ---------------------------------------------------------------------------
  // compress_with_func: Call function on masked elements only
  // ---------------------------------------------------------------------------
  template <typename T, typename F>
  auto compress_with_func(const std::vector<T>& v, const mkfit::Event* ev,
                          F func, const ROOT::RVec<int>& mask) {
    assert(v.size() == mask.size());
    using R = std::decay_t<decltype(func(ev, std::declval<T>()))>;
    ROOT::RVec<R> out;
    const size_t n = v.size();
    out.reserve(n);
    for (size_t i = 0; i < n; ++i) {
      if (mask[i]) out.push_back(func(ev, v[i]));
    }
    return out;
  }

  // ---------------------------------------------------------------------------
  // mask_to_index_vec: Convert boolean mask (1s) to index list
  // e.g., {1, 0, 1, 0, 1} → {0, 2, 4}
  // ---------------------------------------------------------------------------
  inline ROOT::RVec<int> mask_to_index_vec(const ROOT::RVec<int>& mask) {
    ROOT::RVec<int> indices;
    indices.reserve(mask.size() / 16);
    for (size_t i = 0; i < mask.size(); ++i) {
      if (mask[i]) indices.push_back(i);
    }
    return indices;
  }
  // ---------------------------------------------------------------------------
  // neg_mask_to_index_vec: Convert boolean mask (0s) to index list
  // e.g., {1, 0, 1, 0, 1} → {1, 3}
  // ---------------------------------------------------------------------------
  inline ROOT::RVec<int> neg_mask_to_index_vec(const ROOT::RVec<int>& mask) {
    ROOT::RVec<int> indices;
    indices.reserve(mask.size());
    for (size_t i = 0; i < mask.size(); ++i) {
      if (!mask[i]) indices.push_back(i);
    }
    return indices;
  }
}

// Note: both EV_MEMBER and VEC_MEMBER can be data or function, m_id or pT()

#define EV_MAP(EV_MEMBER, VEC_MEMBER) \
  [](const mkfit::Event* ev) { \
    return map_with_member(ev->EV_MEMBER, [](auto& s){ return s.VEC_MEMBER; }); \
  }, { "event" }

#define EV_MAP_FUNC(EV_MEMBER, FUNC) \
  [](const mkfit::Event* ev) { \
    return map_with_func(ev->EV_MEMBER, ev, \
      [](const mkfit::Event* e, const auto& elem) { return e->FUNC(elem); }); \
  }, { "event" }

#define EV_GATHER(EV_MEMBER, VEC_MEMBER, IDX_COLUMN) \
  [](const mkfit::Event* ev, const ROOT::RVec<int> idx_vec) { \
    return gather_with_member(ev->EV_MEMBER, [](auto& s){ return s.VEC_MEMBER; }, idx_vec); \
  }, { "event", IDX_COLUMN }

#define EV_GATHER_FUNC(EV_MEMBER, FUNC, IDX_COLUMN) \
  [](const mkfit::Event* ev, const ROOT::RVec<int> indices) { \
    return gather_with_func(ev->EV_MEMBER, ev, \
      [](const mkfit::Event* e, const auto& elem) { return e->FUNC(elem); }, indices); \
  }, { "event", IDX_COLUMN }

#define EV_COMPRESS(EV_MEMBER, VEC_MEMBER, MASK_COLUMN) \
  [](const mkfit::Event* ev, const ROOT::RVec<int> mask) { \
    return compress_with_member(ev->EV_MEMBER, [](auto& s){ return s.VEC_MEMBER; }, mask); \
  }, { "event", MASK_COLUMN }

#define EV_COMPRESS_FUNC(EV_MEMBER, FUNC, MASK_COLUMN) \
  [](const mkfit::Event* ev, const ROOT::RVec<int> mask) { \
    return compress_with_func(ev->EV_MEMBER, ev, \
      [](const mkfit::Event* e, const auto& elem) { return e->FUNC(elem); }, mask); \
  }, { "event", MASK_COLUMN }


#pragma endregion
//==============================================================================

#endif
