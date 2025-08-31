#include <Selection/CutBase.h>
#include <IsoToolbox/AnalysisContext.h> // 这里需要完整的头文件

namespace PhysicsModules {

CutBase::CutBase(const IsoToolbox::AnalysisContext& context) : m_context(context) {}

}