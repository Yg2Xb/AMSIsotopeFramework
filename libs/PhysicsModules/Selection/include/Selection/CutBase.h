#pragma once

#include <DataModel/StandardizedEvent.h>

namespace IsoToolbox {
    class AnalysisContext;
}

namespace PhysicsModules {

class CutBase {
public:
    CutBase(const IsoToolbox::AnalysisContext& context);
    virtual ~CutBase() = default;

    // The IsPass method now ONLY operates on the StandardizedEvent.
    // Full decoupling from the raw TTree is achieved.
    virtual bool IsPass(const StandardizedEvent& event) const = 0;

protected:
    const IsoToolbox::AnalysisContext& m_context;
};

}