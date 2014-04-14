#ifndef __C3_LINEAR_ALGEBRA_HPP__
#define __C3_LINEAR_ALGEBRA_HPP__

#include <CeedMath/C3Math.h>

#define c3inline inline

#undef c3extern
#define c3extern extern

// Generic Functions
#define C3_DEFINE_GENERIC_SINGLE
#include <CeedMath/C3Math-generic-defs.hpp>
#include <CeedMath/C3Math-generic-protos.h>
#include <CeedMath/C3Math-generic-inline.h>
#undef C3_DEFINE_GENERIC_SINGLE

#define C3_UNDEFINE_GENERIC
#include <CeedMath/C3Math-generic-defs.hpp>
#undef C3_UNDEFINE_GENERIC

#define C3_DEFINE_GENERIC_DOUBLE
#include <CeedMath/C3Math-generic-defs.hpp>
#include <CeedMath/C3Math-generic-protos.h>
#include <CeedMath/C3Math-generic-inline.h>
#undef C3_DEFINE_GENERIC_DOUBLE

#define C3_UNDEFINE_GENERIC
#include <CeedMath/C3Math-generic-defs.hpp>
#undef C3_UNDEFINE_GENERIC

#endif
