#include <Selection/CutBase.h>
#include <IsoToolbox/AnalysisContext.h>

namespace PhysicsModules {

// This is the constructor implementation.
CutBase::CutBase(const IsoToolbox::AnalysisContext& context) : m_context(context) {}

// The IsPass method is pure virtual in the base class, so it has no implementation here.

}