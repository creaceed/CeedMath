#include <CeedMath/C3Math.h>
#include <stdio.h>

#pragma mark globals
const vec3_d kvec3_null_d = {0,0,0};
const vec3_f kvec3_null_f = {0,0,0};
const vec3_d kvec3_xdir_d = {1,0,0};
const vec3_f kvec3_xdir_f = {1,0,0};
const vec3_d kvec3_ydir_d = {0,1,0};
const vec3_f kvec3_ydir_f = {0,1,0};
const vec3_d kvec3_zdir_d = {0,0,1};
const vec3_f kvec3_zdir_f = {0,0,1};

const plane_d kplane_xy_d = {0,0,1,0};
const plane_f kplane_xy_f = {0,0,1,0};
const plane_d kplane_yz_d = {1,0,0,0};
const plane_f kplane_yz_f = {1,0,0,0};
const plane_d kplane_zx_d = {0,1,0,0};
const plane_f kplane_zx_f = {0,1,0,0};


// Generic
#define C3_DEFINE_GENERIC_SINGLE
#include "C3Math-generic-defs.h"
#include "C3Math-generic-extern.h"
#undef C3_DEFINE_GENERIC_SINGLE

#define C3_UNDEFINE_GENERIC
#include "C3Math-generic-defs.h"
#undef C3_UNDEFINE_GENERIC

#define C3_DEFINE_GENERIC_DOUBLE
#include "C3Math-generic-defs.h"
#include "C3Math-generic-extern.h"
#undef C3_DEFINE_GENERIC_DOUBLE

#define C3_UNDEFINE_GENERIC
#include "C3Math-generic-defs.h"
#undef C3_UNDEFINE_GENERIC
