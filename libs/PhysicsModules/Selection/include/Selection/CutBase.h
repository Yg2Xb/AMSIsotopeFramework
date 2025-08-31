#pragma once

// FIX: Include the full definition of the class.
#include <DataModel/AMSDstTreeA.h>

// Forward declaration for AnalysisContext is still appropriate here.
namespace IsoToolbox {
    class AnalysisContext;
}

namespace PhysicsModules {

class CutBase {
public:
    CutBase(const IsoToolbox::AnalysisContext& context);
    virtual ~CutBase() = default;

    // FIX: Removed the incorrect "DataModel::" namespace prefix.
    virtual bool IsPass(AMSDstTreeA* data) const = 0;

protected:
    const IsoToolbox::AnalysisContext& m_context;
};

}