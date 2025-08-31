#pragma once

#include <Selection/CutBase.h> // 必须完整包含父类的头文件

// 子类的头文件同样可以使用前向声明来保持整洁
namespace DataModel {
    class AMSDstTreeA;
}

namespace IsoToolbox {
    class AnalysisContext;
}

namespace PhysicsModules {

class RTICut : public CutBase {
public:
    RTICut(const IsoToolbox::AnalysisContext& context);
    bool IsPass(DataModel::AMSDstTreeA* data) const override;
};

}