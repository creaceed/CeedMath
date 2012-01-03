#ifndef __C3_LINEAR_ALGEBRA_CONVERSIONS_H__
#define __C3_LINEAR_ALGEBRA_CONVERSIONS_H__

// Conversions
c3inline vecn_f vec3_convert_f_from_d(vec3_f res, const vec3_d src)
{
	res[0] = src[0];
	res[1] = src[1];
	res[2] = src[2];
	return res;
}

c3inline matn_f mat4_convert_f_from_d(mat4_f res, const mat4_d src)
{
	int i;
	for(i=0; i<16; i++) res[i] = src[i];
	return res;
}

c3inline frustumn_f frustum_convert_f_from_d(frustum_f res, const frustum_d src)
{
	int i;
	for(i=0; i<8*3; i++) res[i] = src[i];
	return res;
}



#endif