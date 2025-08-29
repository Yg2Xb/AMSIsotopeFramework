#include "Selection/RTICut.h"

namespace Selection {

bool RTICut::pass(const AMSDstTreeA& event) {
    // Logic from old RTICut.h: isGoodRTI(int rtigood)
    return (event.rtigood & 0x7F) == 0;
}

}