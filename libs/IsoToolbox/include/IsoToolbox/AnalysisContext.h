#pragma once
#include "PhysicsConstants.h"

namespace IsoToolbox {

class AnalysisContext {
public:
    AnalysisContext(const std::string& target_particle_name);

    const Isotope& GetTargetIsotope() const { return target_isotope_; }

private:
    const Isotope& target_isotope_;
};

}