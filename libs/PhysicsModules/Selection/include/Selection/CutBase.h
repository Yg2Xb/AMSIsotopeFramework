#pragma once

class AMSDstTreeA;

namespace Selection {

class CutBase {
public:
    virtual ~CutBase(); // <-- REMOVE "= default" from here

    virtual bool pass(const AMSDstTreeA& event) = 0;
};

}