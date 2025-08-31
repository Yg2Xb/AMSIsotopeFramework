#pragma once

// 前向声明 (Forward Declaration)，告诉编译器这些类型存在，但无需知道它们的内部细节
// 这是C++头文件管理的最佳实践
namespace IsoToolbox {
    class AnalysisContext;
}

namespace DataModel {
    class AMSDstTreeA;
}

namespace PhysicsModules {

class CutBase {
public:
    // 构造函数使用 AnalysisContext 的引用，前向声明足够了
    CutBase(const IsoToolbox::AnalysisContext& context);
    virtual ~CutBase() = default;

    // 纯虚函数，定义了所有子类必须实现的接口
    virtual bool IsPass(DataModel::AMSDstTreeA* data) const = 0;

protected:
    const IsoToolbox::AnalysisContext& m_context;
};

}