#pragma once
#include "CutBase.h"

namespace Selection {

class RTICut : public CutBase {
public:
    bool pass(const AMSDstTreeA& event) override;
};

}