#include "RdfVectorSource.h"

#include <TClass.h>
#include <TDataMember.h>

#include <typeinfo>
#include <stdexcept>

namespace mkfit {

  namespace {

    // Helper template to extract type info from TClass members
    template <typename T>
    std::vector<std::string> GetColumnNamesFromClass() {
      std::vector<std::string> names;
      TClass *cls = TClass::GetClass(typeid(T));
      if (!cls) return names;

      TIter next(cls->GetListOfDataMembers());
      TDataMember *member;
      while ((member = dynamic_cast<TDataMember *>(next()))) {
        if (!member->IsaPointer()) {
          names.push_back(member->GetName());
        }
      }
      return names;
    }

    // Helper template to get type name for a specific column
    template <typename T>
    std::string GetTypeNameForColumn(std::string_view columnName) {
      TClass *cls = TClass::GetClass(typeid(T));
      if (!cls) return "";
      TDataMember *member = cls->GetDataMember(columnName.data());
      if (!member) return "";
      return member->GetTypeName();
    }

  }  // anonymous namespace

  // ================================================================
  #pragma region VectorBackedRDataSource Template Implementation
  // ================================================================

  template <typename T>
  class VectorBackedRDataSource final : public ROOT::RDF::RDataSource {
    public:
    using Vec_t = std::vector<T>;

    explicit VectorBackedRDataSource(const Vec_t* vec) : fVec(vec) {}

    std::string GetLabel() override {
        return std::string("VectorBackedRDataSource<") + std::string(typeid(T).name()) + ">";
    }

    // --------------------------------------------------

    std::vector<std::pair<ULong64_t, ULong64_t>>
    GetEntryRanges() override {
      if (!fVec || fNextEntry >= fVec->size()) {
        if (fDebugLevel > 0) {
          printf("VectorBackedRDataSource::GetEntryRanges %s (fNextEntry = %llu, total entries = %zu)\n",
                 fVec ? "end of vector reached" : "vector is null",
                 fNextEntry, fVec ? fVec->size() : 0);
        }
        return {};
      }
      fNextEntry = fVec->size();
      if (fDebugLevel > 0) {
        printf("VectorBackedRDataSource::GetEntryRanges returning range [0, %zu]\n", fVec->size());
      }
      return {{0, fVec->size()}};
    }

    void SetNSlots(unsigned int nSlots) override {
      if (fDebugLevel > 1) {
        printf("VectorBackedRDataSource::SetNSlots called with nSlots = %u\n", nSlots);
      }
      // Chain to the base -- see the note in RdfSources.cc: without it
      // RDataSource::fNSlots stays 0 and RSlotStack::GetSlot() spins forever
      // under implicit MT.
      RDataSource::SetNSlots(nSlots);
      fSlotEntry.resize(nSlots, 0);
    }

    void Initialize() override {
      if (fDebugLevel > 0) {
        printf("VectorBackedRDataSource::Initialize called for %s\n", GetLabel().c_str());
      }
      fNextEntry = 0;
    }
    void InitSlot(unsigned int slot, ULong64_t firstEntry) override {
      if (fDebugLevel > 1) {
        printf("VectorBackedRDataSource::InitSlot called for slot %u, firstEntry %llu\n", slot, firstEntry);
      }
    }
    bool SetEntry(unsigned int slot, ULong64_t entry) override {
      if (fDebugLevel > 2) {
        printf("VectorBackedRDataSource::SetEntry called for slot %u with entry %llu\n", slot, entry);
      }
      fSlotEntry[slot] = entry;
      return true;
    }
    void FinalizeSlot(unsigned int slot) override {
      if (fDebugLevel > 1) {
        printf("VectorBackedRDataSource::FinalizeSlot called for slot %u\n", slot);
      }
    }
    void Finalize() override {
      if (fDebugLevel > 0) {
        printf("VectorBackedRDataSource::Finalize called for %s\n", GetLabel().c_str());
      }
    }

    // --------------------------------------------------
    // Column Introspection via Dictionary
    // --------------------------------------------------

    const std::vector<std::string>& GetColumnNames() const override {
      static std::vector<std::string> names = GetColumnNamesFromClass<T>();
      return names;
    }

    std::string GetTypeName(std::string_view col) const override {
      return GetTypeNameForColumn<T>(col);
    }

    bool HasColumn(std::string_view col) const override {
      auto cls = TClass::GetClass(typeid(T));
      return cls && cls->GetDataMember(col.data());
    }

    // --------------------------------------------------
    // Column Reader
    // --------------------------------------------------

    std::unique_ptr<ROOT::Detail::RDF::RColumnReaderBase>
    GetColumnReaders(unsigned int slot, std::string_view columnName, const std::type_info&) override
    {
      auto cls = TClass::GetClass(typeid(T));
      auto member = cls->GetDataMember(columnName.data());
      const int offset = member->GetOffset();

      auto lambda = [this, slot, offset]() -> void* {
        const auto entry = fSlotEntry[slot];
        if (entry >= fVec->size()) return nullptr;

        return (void*)(
          reinterpret_cast<char*>(
            const_cast<T*>(&(*fVec)[entry])
          ) + offset
        );
      };

      return std::make_unique<LambdaColumnReader>(lambda);
    }

    // Old API - not used but abstract in base
    std::vector<void*> GetColumnReadersImpl(std::string_view, const std::type_info&) override { return {}; }

    private:
    const Vec_t* fVec = nullptr;
    ULong64_t fNextEntry = 0;
    std::vector<ULong64_t> fSlotEntry;
    int fDebugLevel = 1; // 0 = no debug, 1 = basic, 2 = full, 3 = including SetEntry calls
  };

  #pragma endregion
  // ================================================================
  #pragma region RdfVectorSources
  // ================================================================

  ROOT::RDataFrame RdfVectorSources::MakeTrCandMetaDF(const Event &ev) {
    return ROOT::RDataFrame(std::make_unique<VectorBackedRDataSource<TrCandMeta>>(&ev.trCandMetas_));
  }
  ROOT::RDataFrame RdfVectorSources::MakeTrCandStateDF(const Event &ev) {
    return ROOT::RDataFrame(std::make_unique<VectorBackedRDataSource<TrCandState>>(&ev.trCandStates_));
  }
  ROOT::RDataFrame RdfVectorSources::MakeTrHitMatchDF(const Event &ev) {
    return ROOT::RDataFrame(std::make_unique<VectorBackedRDataSource<TrHitMatch>>(&ev.trHitMatches_));
  }

  ROOT::RDataFrame RdfVectorSources::MakeTrackDF(const TrackVec &tvec) {
    return ROOT::RDataFrame(std::make_unique<VectorBackedRDataSource<Track>>(&tvec));
  }
  ROOT::RDataFrame RdfVectorSources::MakeSeedDF(const Event &ev) {
    return MakeTrackDF(*ev.currentSeedTracks_);
  }

  #pragma endregion
}  // end namespace mkfit
