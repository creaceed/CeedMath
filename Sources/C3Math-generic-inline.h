
#pragma mark vec2
// ================== vec2d Functions
vecn_x c3fx(vec2_copy)(vec2_x res, const vec2_x src)
{
	res[0] = src[0];
	res[1] = src[1];
	
	return res;
}
vecn_x c3fx(vec2_add)(vec2_x res, const vec2_x a, const vec2_x b)
{
	res[0] = a[0]+b[0];
	res[1] = a[1]+b[1];
	
	return res;
}
vecn_x c3fx(vec2_subtract)(vec2_x res, const vec2_x a, const vec2_x b)
{
	res[0] = a[0] - b[0];
	res[1] = a[1] - b[1];
	
	return res;
}
vecn_x c3fx(vec2_rotate)(vec2_x res, const vec2_x src, val_x angle)
{
	val_x ca = cos_x(angle), sa = sin_x(angle);
	
	res[0] = src[0]*ca - src[1]*sa;
	res[1] = src[0]*sa + src[1]*ca;
	
	return res;
}
vecn_x c3fx(vec2_lerp)(vec2_x res, const vec2_x a, const vec2_x b, val_x alpha)
{
	val_x oma = 1.0-alpha;
	
	res[0] = a[0]*oma + b[0]*alpha;
	res[1] = a[1]*oma + b[1]*alpha;
	
	return res;
}
val_x c3fx(vec2_square_norm)(const vec2_x a)
{
	return (a[0]*a[0] + a[1]*a[1]);
}
val_x c3fx(vec2_norm)(const vec2_x a)
{
	return sqrt_x(c3fx(vec2_square_norm)(a));
}

val_x c3fx(vec2_distance)(const vec2_x a, const vec2_x b)
{
	vec2_x diff;
	
	return c3fx(vec2_norm)(c3fx(vec2_subtract)(diff,a,b));
}



#pragma mark mat2
// ================== mat2 Functions
matn_x c3fx(mat2_set_rotation)(mat2_x res, val_x a)
{
	val_x ca = cos_x(a), sa = sin_x(a);
	
	res[0] = ca;
	res[1] = -sa;
	res[2] = sa;
	res[3] = ca;
	
	return res;
}
matn_x c3fx(mat2_set_scale)(mat2_x res, val_x sx, val_x sy)
{
	res[0] = sx;
	res[1] = 0;
	res[2] = 0;
	res[3] = sy;
	
	return res;
}
matn_x c3fx(mat2_multiply)(mat2_x res, const mat2_x a, const mat2_x b)
{
	assert(a != res && b != res);
	
	res[0] = a[0]*b[0] + a[1]*b[2];
	res[1] = a[0]*b[1] + a[1]*b[3];
	res[2] = a[2]*b[0] + a[3]*b[2];
	res[3] = a[2]*b[1] + a[3]*b[3];
	
	return res;
}


#pragma mark rect
vecn_x c3fx(rect_point_at_index)(vec2_x res, const rect_x r, int index)
{
	res[0] = r[0] + ((index==1||index==2)?r[2]:0.0);
	res[1] = r[1] + ((index>1)?r[3]:0.0);
	
	return res;
}
vecn_x c3fx(rect_inset)(rect_x res, const rect_x r, val_x dx, val_x dy)
{
	res[0] = r[0] + dx;
	res[1] = r[1] + dy;
	res[2] = r[2] - 2.0*dx;
	res[3] = r[3] - 2.0*dy;
	
	return res;
}

#pragma mark vec3d
// ================== vec3d Functions
vecn_x c3fx(vec3_copy)(vec3_x res, const vec3_x src)
{
	res[0] = src[0];
	res[1] = src[1];
	res[2] = src[2];
	
	return res;
}
vecn_x c3fx(vec3_set)(vec3_x res, val_x x, val_x y, val_x z)
{
	res[0] = x;
	res[1] = y;
	res[2] = z;
	
	return res;
}
vecn_x c3fx(vec3_add)(vec3_x res, const vec3_x a, const vec3_x b)
{
	res[0] = a[0] + b[0];
	res[1] = a[1] + b[1];
	res[2] = a[2] + b[2];
	
	return res;
}
vecn_x c3fx(vec3_subtract)(vec3_x res, const vec3_x a, const vec3_x b)
{
	res[0] = a[0] - b[0];
	res[1] = a[1] - b[1];
	res[2] = a[2] - b[2];
	
	return res;
}
vecn_x c3fx(vec3_negate)(vec3_x res, const vec3_x src)
{
	res[0] = -src[0];
	res[1] = -src[1];
	res[2] = -src[2];
	
	return res;
}
vecn_x c3fx(vec3_scale)(vec3_x res, const vec3_x src, val_x scalar)
{
	res[0] = scalar * src[0];
	res[1] = scalar * src[1];
	res[2] = scalar * src[2];
	
	return res;
}
bool c3fx(vec3_compare)(const vec3_x a, const vec3_x b, val_x tol)
{
	return (fabs_x(b[0]-a[0]) <= tol
			&& fabs_x(b[1]-a[1]) <= tol
			&& fabs_x(b[2]-a[2]) <= tol);
}
bool c3fx(vec3_equal)(const vec3_x a, const vec3_x b)
{
	return c3fx(vec3_compare)(a,b,0.0);
}
bool c3fx(vec3_is_null)(const vec3_x a)
{
	return c3fx(vec3_equal)(a, c3vx(kvec3_null));
}
val_x c3fx(vec3_dot_product)(const vec3_x a, const vec3_x b)
{
	return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}
vecn_x c3fx(vec3_cross_product)(vec3_x res, const vec3_x a, const vec3_x b)
{
	assert(res != a && res != b);
	
	res[0] = a[1] * b[2] - a[2] * b[1];
	res[1] = a[2] * b[0] - a[0] * b[2];
	res[2] = a[0] * b[1] - a[1] * b[0];
	
	return res;
}
val_x c3fx(vec3_mixed_product)(const vec3_x a, const vec3_x b, const vec3_x c)
{
	return (a[0]*b[1]*c[2] + a[1]*b[2]*c[0] + a[2]*b[0]*c[1] - a[2]*b[1]*c[0] - a[0]*b[2]*c[1] - a[1]*b[0]*c[2]);
}
val_x c3fx(vec3_square_norm)(const vec3_x a)
{
	return (a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}
val_x c3fx(vec3_norm)(const vec3_x a)
{
	return sqrt_x(c3fx(vec3_square_norm)(a));
}
val_x c3fx(vec3_distance)(const vec3_x a, const vec3_x b)
{
	vec3_x diff;
	
	return c3fx(vec3_norm)(c3fx(vec3_subtract)(diff,a,b));
}
vecn_x c3fx(vec3_normalize)(vec3_x res, const vec3_x a)
{
	val_x norm = c3fx(vec3_norm)(a);
	
	if(norm > 0.0)
	{
		norm = 1.0 / norm;
		res[0] = a[0] * norm;
		res[1] = a[1] * norm;
		res[2] = a[2] * norm;
	}
	
	return res;
}
vecn_x c3fx(vec3_add_scaled)(vec3_x res, const vec3_x p, const vec3_x q, val_x lambda)
{
	res[0] = p[0] + lambda * q[0];
	res[1] = p[1] + lambda * q[1];
	res[2] = p[2] + lambda * q[2];
	
	return res;
}
vecn_x c3fx(vec3_lerp)(vec3_x res, const vec3_x a, const vec3_x b, val_x alpha)
{
	val_x oma = 1.0-alpha;
	
	res[0] = a[0]*oma + b[0]*alpha;
	res[1] = a[1]*oma + b[1]*alpha;
	res[2] = a[2]*oma + b[2]*alpha;
	
	return res;
}
vecn_x c3fx(vec3_transform)(vec3_x res, const mat4_x mat, const vec3_x v)
{
	vec4_x v4, r4;
	v4[3] = 1.0;
	
	return c3fx(vec3_copy)(res, c3fx(vec4_homogenize)(v4, c3fx(vec4_transform)(r4, mat, c3fx(vec3_copy)(v4, v))));
}
vecn_x c3fx(vec3_project_x)(vec3_x res, const vec3_x a)
{
	res[0] = a[0];
	res[1] = 0.0;
	res[2] = 0.0;
	
	return res;
}
vecn_x c3fx(vec3_project_y)(vec3_x res, const vec3_x a)
{
	res[0] = 0.0;
	res[1] = a[1];
	res[2] = 0.0;
	
	return res;
}
vecn_x c3fx(vec3_project_z)(vec3_x res, const vec3_x a)
{
	res[2] = a[2];
	res[1] = 0.0;
	res[0] = 0.0;
	
	return res;
}
vecn_x c3fx(vec3_project_xy)(vec3_x res, const vec3_x a)
{
	res[0] = a[0];
	res[1] = a[1];
	res[2] = 0.0;
	
	return res;
}
vecn_x c3fx(vec3_project_yz)(vec3_x res, const vec3_x a)
{
	res[2] = a[2];
	res[1] = a[1];
	res[0] = 0.0;
	
	return res;
}
vecn_x c3fx(vec3_project_zx)(vec3_x res, const vec3_x a)
{
	res[0] = a[0];
	res[2] = a[2];
	res[1] = 0.0;
	
	return res;
}

void c3fx(vec3_project_frustum)(vec2_x minmax, const vec3_x dir, const frustum_x frustum)
{
	minmax[0] = MAXFLOAT;
	minmax[1] = -MAXFLOAT;
	
	int i;
	for(i=0; i<8; i++)
	{
		val_x d = c3fx(vec3_dot_product)(frustum+3*i, dir);
		
		minmax[0] = C3MIN(minmax[0], d);
		minmax[1] = C3MAX(minmax[1], d);
	}
}
void c3fx(vec3_project_edge)(vec2_x minmax, const vec3_x dir, const vec3_x at, const vec3_x bt)
{
	minmax[0] = minmax[1] = c3fx(vec3_dot_product)(at, dir);
	val_x d = c3fx(vec3_dot_product)(bt, dir);
	minmax[0] = C3MIN(minmax[0], d);
	minmax[1] = C3MAX(minmax[1], d);
}
void c3fx(vec3_project_triangle)(vec2_x minmax, const vec3_x dir, const vec3_x at, const vec3_x bt, const vec3_x ct)
{
	minmax[0] = minmax[1] = c3fx(vec3_dot_product)(at, dir);
	val_x d = c3fx(vec3_dot_product)(bt, dir);
	minmax[0] = C3MIN(minmax[0], d);
	minmax[1] = C3MAX(minmax[1], d);

	d = c3fx(vec3_dot_product)(ct, dir);
	minmax[0] = C3MIN(minmax[0], d);
	minmax[1] = C3MAX(minmax[1], d);
}

vecn_x c3fx(vec3_center_scale)(vec3_x res, const vec3_x a, const vec3_x center, val_x s)
{
	res[0] = (a[0]-center[0]) * s + center[0];
	res[1] = (a[1]-center[1]) * s + center[1];
	res[2] = (a[2]-center[2]) * s + center[2];
	
	return res;
}
vecn_x c3fx(vec3_reflect)(vec3_x res, const vec3_x a, const normal_x normal)
{
	val_x d = 2.0 * c3fx(vec3_dot_product)(a, normal);
	
	return c3fx(vec3_subtract)(res, a, c3fx(vec3_scale)(res, normal, d));
}
vecn_x c3fx(vec3_refract)(vec3_x res, const vec3_x a, const normal_x normal, val_x index)
{
	val_x d = c3fx(vec3_dot_product)(a, normal), d2 = d*d, r, i = 1.0/index, i2 = i*i;
	
	r = 1.0 - i2*(1.0 - d2);
	if(r < 0.0)
	{
		return c3fx(vec3_reflect)(res, a, normal);
	}
	else
	{
		r = ((d<=0.0)?sqrt_x(r):-sqrt_x(r)) + d*i;
		c3fx(vec3_scale)(res, a, i);
		return c3fx(vec3_add_scaled)(res, res, normal, -r);
	}
}

#pragma mark normal
// ================== normal Functions
bool c3fx(normal_transform)(normal_x res, const mat4_x mat, const normal_x src)
{
	mat4_x mat_it;
	
	if(!c3fx(mat4_invert_transpose)(mat_it, mat)) return false;
	else
	{
		c3fx(normal_transform_it)(res, mat_it, src);
		return true;
	}
}

vecn_x c3fx(normal_transform_it)(normal_x res, const mat4_x mat_it, const normal_x src)
{
	vec4_x v, v2;
	
	c3fx(vec3_copy)(v, src); v[3] = 0.0;
	c3fx(vec4_transform)(v2, mat_it, v);
	return c3fx(vec3_copy)(res, v2);
}

#pragma mark vec4d
// ================== vec4d Functions
vecn_x c3fx(vec4_set)(vec4_x res, val_x x, val_x y, val_x z, val_x w)
{
	res[0] = x;
	res[1] = y;
	res[2] = z;
	res[3] = w;
	
	return res;
}
vecn_x c3fx(vec4_copy)(vec4_x res, const vec4_x src)
{
	res[0] = src[0];
	res[1] = src[1];
	res[2] = src[2];
	res[3] = src[3];
	
	return res;
}

bool c3fx(vec4_compare)(const vec4_x a, const vec4_x b, val_x tol)
{
	return (fabs_x(b[0]-a[0]) <= tol
			&& fabs_x(b[1]-a[1]) <= tol
			&& fabs_x(b[2]-a[2]) <= tol
			&& fabs_x(b[3]-a[3]) <= tol);
}

bool c3fx(vec4_equal)(const vec4_x a, const vec4_x b)
{
	return c3fx(vec4_compare)(a,b,0.0);
}

vecn_x c3fx(vec4_lerp)(vec4_x res, const vec4_x a, const vec4_x b, val_x alpha)
{
	val_x oma = 1.0-alpha;
	
	res[0] = a[0]*oma + b[0]*alpha;
	res[1] = a[1]*oma + b[1]*alpha;
	res[2] = a[2]*oma + b[2]*alpha;
	res[3] = a[3]*oma + b[3]*alpha;
	
	return res;
}

vecn_x c3fx(vec4_homogenize)(vec4_x res, const vec4_x a)
{
	if(fabs_x(a[3]) > 0.0)
	{
		val_x oow = 1.0f/a[3];
		
		res[0] = a[0]*oow;
		res[1] = a[1]*oow;
		res[2] = a[2]*oow;
		res[3] = 1.0;
		
		return res;
	}
	else
	{
		return c3fx(vec4_copy)(res, a);
	}
}
vecn_x c3fx(vec4_transform)(vec4_x res, const mat4_x m, const vec4_x v)
{
	assert(res != v);
	
	res[0] = m[0*4+0]*v[0] + m[0*4+1]*v[1] + m[0*4+2]*v[2] + m[0*4+3]*v[3];
	res[1] = m[1*4+0]*v[0] + m[1*4+1]*v[1] + m[1*4+2]*v[2] + m[1*4+3]*v[3];
	res[2] = m[2*4+0]*v[0] + m[2*4+1]*v[1] + m[2*4+2]*v[2] + m[2*4+3]*v[3];
	res[3] = m[3*4+0]*v[0] + m[3*4+1]*v[1] + m[3*4+2]*v[2] + m[3*4+3]*v[3];
	
	return res;
}
vecn_x c3fx(vec4_transform_delta)(vec4_x res, const mat4_x mat, const vec4_x v)
{
	vec4_x p1, p2, z = {0,0,0,1};
	
	c3fx(vec4_transform)(p1, mat, z);
	c3fx(vec4_transform)(p2, mat, v);
	
	return c3fx(vec3_subtract)(res, p2, p1);
}
vecn_x c3fx(vec4_transpose_transform)(vec4_x res, const mat4_x m, const vec4_x v)
{
	assert(res != v);
	
	res[0] = m[0+4*0]*v[0] + m[0+4*1]*v[1] + m[0+4*2]*v[2] + m[0+4*3]*v[3];
	res[1] = m[1+4*0]*v[0] + m[1+4*1]*v[1] + m[1+4*2]*v[2] + m[1+4*3]*v[3];
	res[2] = m[2+4*0]*v[0] + m[2+4*1]*v[1] + m[2+4*2]*v[2] + m[2+4*3]*v[3];
	res[3] = m[3+4*0]*v[0] + m[3+4*1]*v[1] + m[3+4*2]*v[2] + m[3+4*3]*v[3];
	
	return res;
}


#pragma mark quaternion
// ================== Quaternion Functions
quatn_x c3fx(quat4_set)(quat4_x res, val_x w, val_x x, val_x y, val_x z)
{
	res[0] = w;
	res[1] = x;
	res[2] = y;
	res[3] = z;
	
	return res;
}
quatn_x c3fx(quat4_set_axis_angle)(quat4_x res, const normal_x axis, val_x angle)
{
	val_x sa = sin_x(angle*0.5);
	
	res[0] = cos_x(angle*0.5);
	res[1] = sa*axis[0];
	res[2] = sa*axis[1];
	res[3] = sa*axis[2];
	
	return res;
}

quatn_x c3fx(quat4_set_from_rotation_matrix)(quat4_x res, const mat4_x rotmat)
{
	val_x s = rotmat[MIJ4(0,0)] + rotmat[MIJ4(1,1)] + rotmat[MIJ4(2,2)] + rotmat[MIJ4(3,3)];
	
	res[0] = 1.0; 	res[1] = 0.0; 	res[2] = 0.0; 	res[3] = 0.0; 
	
	if(s > 0.0)
	{
		res[0] = 0.5 * sqrt_x(s);
		res[1] = (rotmat[MIJ4(2,1)] - rotmat[MIJ4(1,2)]) / (4*res[0]);
		res[2] = (rotmat[MIJ4(0,2)] - rotmat[MIJ4(2,0)]) / (4*res[0]);
		res[3] = (rotmat[MIJ4(1,0)] - rotmat[MIJ4(0,1)]) / (4*res[0]);
	}
	
	return res;
}
quatn_x c3fx(quat4_copy)(quat4_x res, const quat4_x src)
{
	return c3fx(vec4_copy)(res, src);
}

quatn_x c3fx(quat4_product)(quat4_x res, const quat4_x q1, const quat4_x q2)
{
	res[0] = q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2] - q1[3]*q2[3];
	res[1] = q1[0]*q2[0] + q1[1]*q2[1] + q1[2]*q2[2] - q1[3]*q2[3];
	res[2] = q1[0]*q2[0] - q1[1]*q2[1] + q1[2]*q2[2] + q1[3]*q2[3];
	res[3] = q1[0]*q2[0] + q1[1]*q2[1] - q1[2]*q2[2] + q1[3]*q2[3];
	
	return res;
}

quatn_x c3fx(quat4_add)(quat4_x res, const quat4_x q1, const quat4_x q2)
{
	res[0] = q1[0]+q2[0];
	res[1] = q1[1]+q2[1];
	res[2] = q1[2]+q2[2];
	res[3] = q1[3]+q2[3];
	
	return res;
}

quatn_x c3fx(quat4_subtract)(quat4_x res, const quat4_x q1, const quat4_x q2)
{
	res[0] = q1[0]-q2[0];
	res[1] = q1[1]-q2[1];
	res[2] = q1[2]-q2[2];
	res[3] = q1[3]-q2[3];
	
	return res;	
}

quatn_x c3fx(quat4_scale)(quat4_x res, const quat4_x q, val_x t)
{
	res[0] = q[0]*t;
	res[1] = q[1]*t;
	res[2] = q[2]*t;
	res[3] = q[3]*t;
	
	return res;	
}
val_x    c3fx(quat4_norm)(const quat4_x q)
{
	return sqrt_x(C3SQR(q[0])+C3SQR(q[1])+C3SQR(q[2])+C3SQR(q[3]));
}

quatn_x c3fx(quat4_normalize)(quat4_x res, const quat4_x q)
{
	val_x n = c3fx(quat4_norm)(q);
	
	if(n == 0.0) n=1.0;
	n = 1.0/n;
	
	return c3fx(quat4_scale)(res, q, n);
}

quatn_x c3fx(quat4_slerp)(quat4_x res, const quat4_x q1, const quat4_x q2, val_x t)
{
	quat4_x sum, dif, q2p;
	val_x dp, s, a1, a2;
	
	c3fx(quat4_copy)(res, q1);
	c3fx(quat4_copy)(q2p, q2);
	c3fx(quat4_add)(sum, q1, q2);
	c3fx(quat4_subtract)(dif, q1, q2);
	
	if(c3fx(quat4_norm)(dif) >= c3fx(quat4_norm)(sum)) c3fx(quat4_scale)(q2p, q2, -1.0);
	
	dp = q1[0]*q2p[0] + q1[1]*q2p[1] + q1[2]*q2p[2] + q1[3]*q2p[3];
	
	if(fabs_x(dp) <= 1.0)
	{
		dp = acos_x(dp);
		s = sin_x(dp);
		
		if(fabs_x(s)>0.0)
		{
			a1 = sin_x((1.0-t)*dp) / s;
			a2 = sin_x(t*dp) / s;
			
			res[0] = q1[0]*a1 + q2[0]*a2;
			res[1] = q1[1]*a1 + q2[1]*a2;
			res[2] = q1[2]*a1 + q2[2]*a2;
			res[3] = q1[3]*a1 + q2[3]*a2;
		}
	}
	
	return res;
}

#pragma mark plane
// ================== Plane Functions

vecn_x c3fx(plane_copy)(plane_x res, const plane_x src)
{
	res[0] = src[0];
	res[1] = src[1];
	res[2] = src[2];
	res[3] = src[3];
	
	return res;
}

vecn_x c3fx(plane_set_normal_point)(plane_x res, const normal_x normal, const vec3_x p)
{
	c3fx(vec3_copy)(res, normal);
	res[3] = - c3fx(vec3_dot_product)(res, p);
	
	return res;
}

vecn_x c3fx(plane_set_points)(plane_x res, const vec3_x a, const vec3_x b, const vec3_x c)
{
	vec3_x ab, bc, n;
	
	c3fx(vec3_subtract)(ab, b, a);
	c3fx(vec3_subtract)(bc, c, b);
	c3fx(vec3_cross_product)(n, ab, bc);
	
	return c3fx(plane_set_normal_point)(res, n, a);
}
vecn_x c3fx(plane_transform_it)(plane_x res, const mat4_x mat_it, const plane_x src)
{
	return c3fx(vec4_transform)(res, mat_it, src);
}

bool c3fx(plane_transform)(plane_x res, const mat4_x mat, const plane_x src)
{
	mat4_x mat_it;
	
	if(c3fx(mat4_invert_transpose)(mat_it, mat))
	{
		c3fx(plane_transform_it)(res, mat_it, src);
		return true;
	}
	
	return false;
}


#pragma mark mat4d
// ================== Matrix 4 Functions
matn_x c3fx(mat4_copy)(mat4_x res, const mat4_x src)
{
	memcpy(res, src, 16 * sizeof(res[0]));
	
	return res;
}
matn_x c3fx(mat4_multiply)(mat4_x res, const mat4_x a, const mat4_x b)
{
	assert(res != a && res != b);
	
	int i,j;
	
	for(i=0; i<4; i++)
	{
		for(j=0; j<4; j++)
		{
			res[MIJ4(i,j)] = a[MIJ4(i,0)]*b[MIJ4(0,j)] + a[MIJ4(i,1)]*b[MIJ4(1,j)] + a[MIJ4(i,2)]*b[MIJ4(2,j)] + a[MIJ4(i,3)]*b[MIJ4(3,j)];
		}
	}
	
	return res;
}


#pragma mark ray
// ================== ray Functions
rayn_x c3fx(ray_set)(ray_x res, const vec3_x origin, const vec3_x dir)
{
	res[0] = origin[0]; 
	res[1] = origin[1]; 
	res[2] = origin[2];
	
	res[3] = dir[0]; 
	res[4] = dir[1]; 
	res[5] = dir[2];
	
	return res;
}
rayn_x c3fx(ray_set_xaxis)(ray_x res)
{
	res[0] = 0.0; 
	res[1] = 0.0; 
	res[2] = 0.0;
	
	res[3] = 1.0; 
	res[4] = 0.0; 
	res[5] = 0.0;
	
	return res;
}
rayn_x c3fx(ray_set_yaxis)(ray_x res)
{
	res[0] = 0.0; 
	res[1] = 0.0; 
	res[2] = 0.0;
	
	res[3] = 0.0; 
	res[4] = 1.0; 
	res[5] = 0.0;
	
	return res;
}
rayn_x c3fx(ray_set_zaxis)(ray_x res)
{
	res[0] = 0.0; 
	res[1] = 0.0; 
	res[2] = 0.0;
	
	res[3] = 0.0; 
	res[4] = 0.0; 
	res[5] = 1.0;
	
	return res;
}
rayn_x c3fx(ray_copy)(ray_x res, ray_x src)
{
	memcpy(res, src, sizeof(ray_d));
	
	return res;
}
vecn_x c3fx(ray_get_origin)(vec3_x res, const ray_x ray)
{
	res[0] = ray[0];
	res[1] = ray[1];
	res[2] = ray[2];
	
	return res;
}
vecn_x c3fx(ray_get_direction)(vec3_x res, const ray_x ray)
{
	res[0] = ray[3];
	res[1] = ray[4];
	res[2] = ray[5];
	
	return res;
}



