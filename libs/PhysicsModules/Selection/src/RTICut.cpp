#include <Selection/RTICut.h>
#include <IsoToolbox/AnalysisContext.h>

namespace PhysicsModules {

RTICut::RTICut(const IsoToolbox::AnalysisContext& context) : CutBase(context) {}

bool RTICut::IsPass(const StandardizedEvent& event) const {
    // The implementation now uses the clean variables from the StandardizedEvent struct.
    
    // For MC, we typically bypass RTI cuts
    if (event.is_mc) {
        return true;
    }

    // Example using a clean, standardized variable:
    if (event.is_bad_run) {
        return false;
    }
    
    // Additional cut logic would go here, using only members of 'event'.
    // Example: if(event.zenith > 40) return false;
    // (Note: 'zenith' would need to be added to StandardizedEvent first)
    
    return true; // Placeholder
}

}