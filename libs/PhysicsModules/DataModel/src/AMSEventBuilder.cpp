#include "DataModel/AMSEventBuilder.h"
#include "DataModel/AMSDstTreeA.h" // Raw TTree structure is now a private implementation detail
#include "IsoToolbox/Logger.h"
#include "IsoToolbox/Tools.h"
#include <algorithm> // For std::copy

namespace PhysicsModules {

// PImpl (Pointer to Implementation) Idiom
// This struct contains the raw TTree data object. By defining it here in the .cpp file,
// we completely hide the AMSDstTreeA class from any code that just includes AMSEventBuilder.h.
class AMSEventBuilder::RawData {
public:
    std::unique_ptr<AMSDstTreeA> event;
    RawData() : event(std::make_unique<AMSDstTreeA>()) {}
};

AMSEventBuilder::AMSEventBuilder(TTree* tree, const IsoToolbox::AnalysisContext* context)
    : m_rawData(std::make_unique<RawData>()) {
    
    if (!tree) {
        throw std::runtime_error("AMSEventBuilder: TTree pointer is null.");
    }
    m_rawData->event->Init(tree);
    MapConfigToAccessors(context);
    IsoToolbox::Logger::Info("AMSEventBuilder initialized successfully.");
}

AMSEventBuilder::~AMSEventBuilder() = default;

void AMSEventBuilder::MapConfigToAccessors(const IsoToolbox::AnalysisContext* context) {
    const auto& chain = context->GetAnalysisChain();

    if (chain.rigidity_version == "GBL_v6_InnerSpan") {
        m_rigidityAccessor      = [](const RawData* rd) { return rd->event->tk_rigidity1[1][2][1]; };
        m_rigidityErrorAccessor = [](const RawData* rd) { return rd->event->tk_erigidity1[1][1]; };
        IsoToolbox::Logger::Info("  -> Rigidity mapped to: GBL_v6_InnerSpan (tk_rigidity1[1][2][1])");
    } else {
        m_rigidityAccessor      = [](const RawData* rd) { return rd->event->tk_rigidity1[1][2][1]; };
        m_rigidityErrorAccessor = [](const RawData* rd) { return rd->event->tk_erigidity1[1][1]; };
        IsoToolbox::Logger::Warn("  -> Rigidity version '{}' not recognized, using default GBL_v6_InnerSpan.", chain.rigidity_version);
    }
    
    if (chain.beta_version == "TOF_BetaH_v1") {
        m_betaTofAccessor = [](const RawData* rd) { return rd->event->tof_betah; };
        IsoToolbox::Logger::Info("  -> TOF Beta mapped to: TOF_BetaH_v1 (tof_betah)");
    } else {
        m_betaTofAccessor = [](const RawData* rd) { return rd->event->tof_betah; };
        IsoToolbox::Logger::Warn("  -> TOF Beta version '{}' not recognized, using default TOF_BetaH_v1.", chain.beta_version);
    }
}

void AMSEventBuilder::Fill(StandardizedEvent& event) const {
    const auto* raw = m_rawData->event.get();

    // --- General Event Information ---
    event.run = raw->run;
    event.event = raw->event;
    event.is_mc = !raw->isreal;
    event.is_bad_run = raw->isbadrun;
    event.cut_status = raw->cutStatus;

    // --- Standardized & Configured Physics Quantities ---
    event.rigidity = m_rigidityAccessor(m_rawData.get());
    event.rigidity_error = m_rigidityErrorAccessor(m_rawData.get());
    event.beta_tof = m_betaTofAccessor(m_rawData.get());
    
    // --- Directly Mapped Common Raw Branches ---
    event.rigidity_inner = raw->tk_rigidity1[1][2][1];
    event.rigidity_l1 = raw->tk_rigidity1[1][2][2];
    event.beta_rich_naf = raw->rich_NaF ? raw->rich_beta[0] : -1.0f;
    event.beta_rich_agl = !raw->rich_NaF ? raw->rich_beta[0] : -1.0f;
    event.ek_tof = raw->TOFEk;
    event.ek_rich_naf = raw->NaFEk;
    event.ek_rich_agl = raw->AGLEk;
    event.cutoff = raw->cutOffRig;
    
    event.charge_l1_unbiased = raw->tk_exqln[0][0][2];
    event.charge_l1_normal = raw->tk_qln[0][0][2];
    event.charge_l2_normal = raw->tk_qln[0][1][2];
    event.charge_inner = raw->tk_qin[0][2];
    event.charge_tof_upper = 0.5f * (raw->tof_ql[0] + raw->tof_ql[1]);
    event.charge_tof_lower = 0.5f * (raw->tof_ql[2] + raw->tof_ql[3]);
    event.charge_rich = raw->rich_q[0];

    event.rich_pmt = raw->rich_pmt;
    event.rich_hit = raw->rich_hit;
    event.rich_used_hits = raw->rich_used;
    std::copy(std::begin(raw->rich_npe), std::end(raw->rich_npe), std::begin(event.rich_npe));
    event.rich_is_naf = raw->rich_NaF;
    event.rich_beta_consistency = raw->rich_BetaConsistency;
    
    event.ntrack = raw->ntrack;
    std::copy(std::begin(raw->tk_qhit), std::end(raw->tk_qhit), std::begin(event.tk_q_hit));
    std::copy(std::begin(raw->tk_qrms), std::end(raw->tk_qrms), std::begin(event.tk_q_rms));
    event.tk_inner_q_rms = raw->tk_qrmn[0][2];

    if (event.is_mc) {
        event.generated_rigidity = raw->mmom / raw->mch;
        event.generated_ek = IsoToolbox::Tools::rigidityToEkPerNucleon(event.generated_rigidity, raw->mpar/1000, raw->mpar%1000);
        event.particle_id_L1 = raw->mtrpar[0];
        event.particle_id_L2 = raw->mtrpar[1];
    } else {
        event.generated_rigidity = 0.0;
        event.generated_ek = 0.0;
        event.particle_id_L1 = 0;
        event.particle_id_L2 = 0;
    }
}

} // namespace PhysicsModules