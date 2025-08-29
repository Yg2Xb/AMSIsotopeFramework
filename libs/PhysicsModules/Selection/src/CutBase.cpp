#include "Selection/CutBase.h"

// By providing an implementation in a .cpp file, we give the linker a
// single, unambiguous location for the class's virtual table,
// which resolves the "undefined reference to vtable" error.
Selection::CutBase::~CutBase() = default;