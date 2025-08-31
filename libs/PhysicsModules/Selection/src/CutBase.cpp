#pragma once

// Forward declaration is efficient here
namespace IsoToolbox {
    class AnalysisContext;
}

// Forward declaration
namespace DataModel {
    class AMSDstTreeA;
}

namespace PhysicsModules {

class CutBase {
public:
    // FIX: Added the necessary forward declaration for AnalysisContext above
    CutBase(const IsoToolbox::AnalysisContext& context);
    virtual ~CutBase() = default;

    virtual bool IsPass(DataModel::AMSDstTreeA* data) const = 0;

protected:
    const IsoToolbox::AnalysisContext& m_context;
};

}