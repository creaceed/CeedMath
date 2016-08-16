#ifndef __C3_LINEAR_ALGEBRA_H__
#define __C3_LINEAR_ALGEBRA_H__

#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <assert.h>


#ifdef __cplusplus
#define c3extern extern "C"
#else
#define c3extern extern
#endif



#ifdef __APPLE__
  #include <TargetConditionals.h>
  #if TARGET_CPU_ARM == 1
	#define c3inline static __inline__ __attribute__((always_inline))
//#define c3inline static __inline__ __attribute__((always_inline, used))
  #else
    // this is the same as NS_INLINE (but still generates warnings)
    #define c3inline static __inline__ __attribute__((always_inline, used))
    //#define c3inline __attribute__((always_inline, used))
    //#define c3inline __attribute__((always_inline))
  #endif
#else
  #define c3inline static inline
#endif

typedef double val_d;
typedef float val_f;

typedef val_d vec2_d[2];
typedef val_f vec2_f[2];
typedef val_d vec3_d[3];
typedef val_f vec3_f[3];
typedef val_d vec4_d[4];
typedef val_f vec4_f[4];
typedef val_d *vecn_d; 	// return value
typedef val_f *vecn_f;

typedef val_d quat4_d[4];
typedef val_f quat4_f[4];
typedef val_d *quatn_d;
typedef val_f *quatn_f;

typedef val_d mat2_d[4];
typedef val_f mat2_f[4];
typedef val_d mat3_d[9];
typedef val_f mat3_f[9];
typedef val_d mat4_d[16];
typedef val_f mat4_f[16];
typedef val_d *matn_d;		// return value
typedef val_f *matn_f;

typedef vec3_d normal_d;	// normalized variant of vec3
typedef vec3_f normal_f;
typedef vec4_d plane_d;	// plane equation A.x + B.y + C.z + D = 0
typedef vec4_f plane_f;

typedef val_d ray_d[6];		// ray in 3d, (ox, oy, oz), (dx, dy, dz), direction not necessarily normalized
typedef val_f ray_f[6];
typedef val_d *rayn_d;		// return value
typedef val_f *rayn_f;

typedef val_d rect_d[4];		// 2d rectangle, (x,y,w,h)
typedef val_f rect_f[4];

typedef val_d box_d[6];	// first 3 for lower corner, last 3 for upper one
typedef val_f box_f[6];
typedef val_d *boxn_d;		// return value
typedef val_f *boxn_f;

typedef val_d frustum_d[8*3];	// first 4 are near, then far, in that order: bottom-left, bottom-right, top-right, top-left
typedef val_f frustum_f[8*3];
typedef val_d *frustumn_d;	// return value
typedef val_f *frustumn_f;

// C3 self contained structs, and light weight conversions
typedef union _vec3d_s {
	struct { val_d x, y, z; };
	vec3_d comps;
} CSVector;
#define CSVectorFromVec3d(vt) (*(CSVector*)vt)

typedef union _vec3f_s {
	struct { val_f x, y, z; };
	vec3_f comps;
} CSVectorf;
#define CSVectorFromVec3d(vt) (*(CSVector*)vt)

typedef union _rayd_s {
	struct {
		struct { val_d x,y,z; } origin;
		struct { val_d x,y,z; } direction;
	};
	ray_d comps;
} CSRay;
#define CSRayFromRay(vt) (*(CSRay*)vt)

typedef union _plane_s {
	struct {
		struct { val_d a,b,c,d; };
	};
	vec4_d comps;
} CSPlane;
#define CSPlaneFromPlane(vt) (*(CSPlane*)vt)

typedef union _quaternion_s {
	struct {
		struct { val_d w,x,y,z; };
	};
	vec4_d comps;
} CSQuaternion;

typedef union _frustum_s {
	struct {
		CSVector near[4]; // bottom-left, bottom-right, top-right, top-left
		CSVector far[4]; // idem
	};
	frustum_d comps;
} CSFrustum;

#define CSQuaternionFromQuat(vt) (*(CSQuaternion*)vt)

// Quickies
#define MIJ4(i,j) ((i)*4+(j))
#define MIJ3(i,j) ((i)*3+(j))
#define C3MIN(a,b) ((a)<(b)?(a):(b))
#define C3MAX(a,b) ((a)>(b)?(a):(b))
#define C3SQR(a) ((a)*(a))

#define C3_RAY_ORIGIN(ray) (ray)
#define C3_RAY_DIRECTION(ray) (&((ray)[3]))

// Some math macros
#define C3_LERP(a,b,r) ((r)*((b)-(a)) + (a))
#define C3_CLAMP(x, min, max) (((x)<(min))?(min):(((x)>(max))?(max):(x)))
#define C3_DEG2RAD(xxx) ((xxx)*M_PI/180.0)
#define C3_RAD2DEG(xxx) ((xxx)/M_PI*180.0)


// globals
extern const vec3_d kvec3_null_d;
extern const vec3_f kvec3_null_f;
extern const vec3_d kvec3_xdir_d;
extern const vec3_f kvec3_xdir_f;
extern const vec3_d kvec3_ydir_d;
extern const vec3_f kvec3_ydir_f;
extern const vec3_d kvec3_zdir_d;
extern const vec3_f kvec3_zdir_f;

extern const plane_d kplane_xy_d;
extern const plane_f kplane_xy_f;
extern const plane_d kplane_yz_d;
extern const plane_f kplane_yz_f;
extern const plane_d kplane_zx_d;
extern const plane_f kplane_zx_f;

#include <CeedMath/C3Math-conversions.h>

// Generic Functions
#define C3_DEFINE_GENERIC_SINGLE
#include <CeedMath/C3Math-generic-defs.h>
#include <CeedMath/C3Math-generic-protos.h>
#include <CeedMath/C3Math-generic-inline.h>
#undef C3_DEFINE_GENERIC_SINGLE

#define C3_UNDEFINE_GENERIC
#include <CeedMath/C3Math-generic-defs.h>
#undef C3_UNDEFINE_GENERIC

#define C3_DEFINE_GENERIC_DOUBLE
#include <CeedMath/C3Math-generic-defs.h>
#include <CeedMath/C3Math-generic-protos.h>
#include <CeedMath/C3Math-generic-inline.h>
#undef C3_DEFINE_GENERIC_DOUBLE

#define C3_UNDEFINE_GENERIC
#include <CeedMath/C3Math-generic-defs.h>
#undef C3_UNDEFINE_GENERIC

#endif
