#ifndef AMS_EVENT_BUILDER_H
#define AMS_EVENT_BUILDER_H

#include "DataModel/StandardizedEvent.h"
#include "IsoToolbox/AnalysisContext.h"
#include <TTree.h>
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <functional>

// Forward declaration of the internal raw data structure
class AMSDstTreeA;

namespace PhysicsModules {

// The core engine of the Data Access Layer. Its sole responsibility is to
// translate data from the raw AMSDstTreeA into the comprehensive StandardizedEvent format.
class AMSEventBuilder {
public:
    AMSEventBuilder(TTree* tree, const IsoToolbox::AnalysisContext* context);
    ~AMSEventBuilder();

    // Fills the provided StandardizedEvent struct with data from the current TTree entry.
    void Fill(StandardizedEvent& event) const;

private:
    void MapConfigToAccessors(const IsoToolbox::AnalysisContext* context);
    
    // PImpl idiom: Pointer to a private, forward-declared struct that holds
    // the raw TTree data. This completely hides AMSDstTreeA from any file
    // that includes AMSEventBuilder.h.
    class RawData;
    std::unique_ptr<RawData> m_rawData;
    
    // Configurable Accessors
    std::function<float(const RawData*)> m_rigidityAccessor;
    std::function<float(const RawData*)> m_rigidityErrorAccessor;
    std::function<float(const RawData*)> m_betaTofAccessor;
};

} // namespace PhysicsModules

#endif // AMS_EVENT_BUILDER_H