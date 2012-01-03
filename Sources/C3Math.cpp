#include "C3Math.hpp"
#include <stdio.h>

// Generic
#define C3_DEFINE_GENERIC_SINGLE
#include "C3Math-generic-defs.hpp"
#include "C3Math-generic-extern.h"
#undef C3_DEFINE_GENERIC_SINGLE

#define C3_UNDEFINE_GENERIC
#include "C3Math-generic-defs.hpp"
#undef C3_UNDEFINE_GENERIC

#define C3_DEFINE_GENERIC_DOUBLE
#include "C3Math-generic-defs.hpp"
#include "C3Math-generic-extern.h"
#undef C3_DEFINE_GENERIC_DOUBLE

#define C3_UNDEFINE_GENERIC
#include "C3Math-generic-defs.hpp"
#undef C3_UNDEFINE_GENERIC
