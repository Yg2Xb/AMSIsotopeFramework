#pragma once

#include <Selection/CutBase.h> // This now brings in the correct definition

namespace IsoToolbox {
    class AnalysisContext;
}

namespace PhysicsModules {

class RTICut : public CutBase {
public:
    RTICut(const IsoToolbox::AnalysisContext& context);
    
    // FIX: Removed the incorrect "DataModel::" namespace prefix to match the base class.
    bool IsPass(AMSDstTreeA* data) const override;
};

}