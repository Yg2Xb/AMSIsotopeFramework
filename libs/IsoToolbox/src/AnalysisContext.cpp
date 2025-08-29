#include "IsoToolbox/AnalysisContext.h"

namespace IsoToolbox {

AnalysisContext::AnalysisContext(const std::string& target_particle_name)
    : target_isotope_(PhysicsConstants::GetIsotope(target_particle_name))
{}

}