#include <cmath>

#ifdef C3_DEFINE_GENERIC_SINGLE
#define c3vx(name) name ## _f
#endif
#ifdef C3_DEFINE_GENERIC_DOUBLE
#define c3vx(name) name ## _d
#endif

#if defined (C3_DEFINE_GENERIC_SINGLE) || defined (C3_DEFINE_GENERIC_DOUBLE)

#define c3fx(name) name	
#define sin_x sin
#define cos_x cos
#define tan_x tan
#define asin_x asin
#define acos_x acos
#define atan2_x atan2
#define sqrt_x sqrt
#define fabs_x fabs
#define floor_x floor

#define val_x  c3vx(val)
#define vecn_x c3vx(vecn)
#define vec2_x c3vx(vec2)
#define vec3_x c3vx(vec3)
#define vec4_x c3vx(vec4)
#define matn_x c3vx(matn)
#define mat2_x c3vx(mat2)
#define mat3_x c3vx(mat3)
#define mat4_x c3vx(mat4)
#define quatn_x c3vx(quatn)
#define quat4_x c3vx(quat4)
#define normal_x c3vx(normal)
#define plane_x c3vx(plane)
#define rayn_x c3vx(rayn)
#define ray_x c3vx(ray)
#define rect_x c3vx(rect)
#define box_x c3vx(box)
#define boxn_x c3vx(boxn)
#define frustum_x c3vx(frustum)
#define frustumn_x c3vx(frustumn)

#endif

// Undefined
#ifdef C3_UNDEFINE_GENERIC
#undef sin_x
#undef cos_x
#undef tan_x
#undef asin_x
#undef acos_x
#undef atan2_x
#undef sqrt_x
#undef fabs_x
#undef floor_x

#undef c3fx
#undef c3vx

#undef val_x
#undef vecn_x
#undef vec2_x
#undef vec3_x
#undef vec4_x
#undef matn_x
#undef mat2_x
#undef mat3_x
#undef mat4_x
#undef quatn_x
#undef quat4_x
#undef normal_x
#undef plane_x
#undef rayn_x
#undef ray_x
#undef rect_x
#undef box_x
#undef boxn_x
#undef frustum_x
#undef frustumn_x
#endif
