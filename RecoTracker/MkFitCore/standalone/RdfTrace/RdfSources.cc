#include "RdfSources.h"

#include <TClass.h>
#include <TDataMember.h>

#include <typeinfo>
#include <stdexcept>

namespace mkfit {

  // ================================================================
  #pragma region EventSource
  // ================================================================

  /**
   * @brief A custom RDataSource that wraps a single Event or a collection of Events.
   *
   * Design: One RDF entry = one mkfit::Event.
   * The data source exposes a single column "event" of type const mkfit::Event*.
   * All other data (seeds, hits, etc.) is accessed via lambda projections on this pointer.
   */
  class EventSource : public ROOT::RDF::RDataSource {
  public:
    explicit EventSource(std::vector<const Event*>& events) : fEvents(events) {
      if (!fEvents.empty()) {
        fNextEntry = 0;
      }
    }

    ~EventSource() override = default;

    std::string GetLabel() override {
      return "EventSource";
    }

    // --------------------------------------------------

    std::vector<std::pair<ULong64_t, ULong64_t>>
    GetEntryRanges() override {
      ULong64_t ev_size = fEvents.size();

      if (fEvents.empty() || fNextEntry >= ev_size) {
        if (fDebugLevel > 0) {
          printf("EventSource::GetEntryRanges called ... vector is null or fNextEntry out of bounds\n");
        }
        return {};
      }
      std::vector<std::pair<ULong64_t, ULong64_t>> out;
      ULong64_t ns = fSlotEntry.size();
      ULong64_t nentries = ev_size - fNextEntry;
      ULong64_t atom = std::max(nentries / ns, 1ull);
      while (fNextEntry < ev_size) {
        ULong64_t beg = fNextEntry;
        ULong64_t end = std::min(beg + atom, ev_size);
        if (fDebugLevel > 0) {
          printf("EventSource::GetEntryRanges appending range [%llu, %llu]\n", beg, end);
        }
        out.push_back( { beg, end } );
        fNextEntry = end;
      }
      return out;
    }

    void SetNSlots(unsigned int nSlots) override {
      if (fDebugLevel > 1) {
        printf("EventSource::SetNSlots called with nSlots = %u\n", nSlots);
      }
      // MUST chain to the base: RDataSource::ProcessMT() sizes its RSlotStack
      // from RDataSource::fNSlots, which only the base setter assigns. Leaving
      // it at 0 gives an empty slot stack, and RSlotStack::GetSlot() then spins
      // forever in its `while (true)` -- an instant, silent deadlock under
      // implicit MT.
      RDataSource::SetNSlots(nSlots);
      fSlotEntry.resize(nSlots, nullptr);
    }

    void Initialize() override {
      fNextEntry = 0;
      if (fDebugLevel > 0) {
        printf("EventSource::Initialize called for %s\n", GetLabel().c_str());
      }
    }

    void InitSlot(unsigned int slot, ULong64_t firstEntry) override {
      if (fDebugLevel > 1) {
        printf("EventSource::InitSlot called for slot %u, firstEntry %llu\n", slot, firstEntry);
      }
    }

    bool SetEntry(unsigned int slot, ULong64_t entry) override {
      if (fDebugLevel > 2) {
        printf("EventSource::SetEntry called for slot %u with entry %llu\n", slot, entry);
      }
      if (entry >= fEvents.size()) return false;
      fSlotEntry[slot] = fEvents[entry];
      return true;
    }

    // NOTE (ROOT dev-1, seen 2026-09-04 with fDebugLevel = 3): under implicit MT
    // this is called TWICE per slot, and in the serial case there is a spurious
    // InitSlot/FinalizeSlot pair around the final, empty GetEntryRanges().
    // Harmless while these are no-ops -- but do NOT put real teardown here
    // (releasing an REve object, say) without re-checking against a newer ROOT.
    void FinalizeSlot(unsigned int slot) override {
      if (fDebugLevel > 1) {
          printf("EventSource::FinalizeSlot called for slot %u\n", slot);
      }
    }

    void Finalize() override {if (fDebugLevel > 0) {
        printf("EventSource::Finalize called for %s\n", GetLabel().c_str());
      }
    }

    // --------------------------------------------------
    // Column Introspection -- trivial
    // --------------------------------------------------

    const std::vector<std::string>& GetColumnNames() const override {
      static std::vector<std::string> names = { "event" };
      return names;
    }

    std::string GetTypeName(std::string_view col) const override {
      return "const mkfit::Event*";
    }

    bool HasColumn(std::string_view col) const override {
      return (col == "event");
    }

    // --------------------------------------------------
    // Column Reader
    // --------------------------------------------------

    std::unique_ptr<ROOT::Detail::RDF::RColumnReaderBase>
    GetColumnReaders(unsigned int slot, std::string_view columnName, const std::type_info&) override {
      if (columnName != "event")
        return nullptr;

      auto lambda = [this, slot]() -> void* {
        return (void*) &fSlotEntry[slot];
      };

      return std::make_unique<LambdaColumnReader>(lambda);
    }

    // Old API - not used but abstract in base
    std::vector<void*> GetColumnReadersImpl(std::string_view, const std::type_info&) override { return {}; }

  private:
    std::vector<const Event*> &fEvents;
    std::vector<const Event*>  fSlotEntry;
    ULong64_t fNextEntry = 0;
    int fDebugLevel = 0; // 0 = no debug, 1 = basic, 2 = full, 3 = including SetEntry calls
  };

  #pragma endregion
  // ================================================================
  #pragma region RdfSources
  // ================================================================

  ROOT::RDataFrame RdfSources::MakeEventDF(std::vector<const Event*>& events) {
    if (events.empty()) return ROOT::RDataFrame(0);
    return ROOT::RDataFrame(std::make_unique<EventSource>(events));
  }

  #pragma endregion
}  // end namespace mkfit
