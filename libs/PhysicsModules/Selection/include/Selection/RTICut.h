#pragma once

#include <Selection/CutBase.h> 

namespace IsoToolbox {
    class AnalysisContext;
}

namespace PhysicsModules {

class RTICut : public CutBase {
public:
    RTICut(const IsoToolbox::AnalysisContext& context);
    
    // --- MODIFIED ---
    // The IsPass override now correctly matches the new base class signature.
    bool IsPass(const StandardizedEvent& event) const override;
};

}
