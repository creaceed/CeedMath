
#pragma mark vec3
vecn_x c3fx(vec3_rotate_centered_axis)(vec3_x res, const vec3_x src, val_x a, const normal_x axis, const vec3_x center)
{
	mat4_x mat;
	vec3_x tmp;
	
	c3fx(mat4_set_rotation_axis_angle)(mat, axis, a);
	c3fx(vec3_subtract)(tmp, src, center);
	c3fx(vec3_transform)(res, mat, tmp);
	
	return c3fx(vec3_add)(res, res, center);
}
val_x c3fx(vec3_angle)(const vec3_x a, const vec3_x b)
{
	val_x n1, n2;
	
	n1 = c3fx(vec3_norm)(a);
	n2 = c3fx(vec3_norm)(b);
	
	if(n1 > 0.0 && n2 > 0.0)
	{
		val_x r = c3fx(vec3_dot_product)(a,b)/(n1*n2);
		
		if(r>1.0) r = 1.0;
		else if(r<-1.0) r = -1.0;
		
		return acos_x(r);
	}
	
	else return 0.0;
}
val_x c3fx(vec3_oriented_angle)(const vec3_x a, const vec3_x b, const vec3_x ref)
{
	val_x n;
	
	n = c3fx(vec3_square_norm)(a) * c3fx(vec3_square_norm)(b);
	
	if(n > 0.0)
	{
		vec3_x cr;
		val_x r = c3fx(vec3_dot_product)(a, b) / sqrt_x(n);
		val_x an;
		
		if(r>1.0) r=1.0;
		else if(r<-1.0) r=-1.0;
		
		an = acos_x(r);
		
		c3fx(vec3_cross_product)(cr, a, b);
		if(c3fx(vec3_dot_product)(cr, ref) >= 0.0) return -an;
		else return an;
	}
	
	else return 0.0;
}

#pragma mark vec4
bool c3fx(vec4_solve_system)(vec4_x res, const mat4_x mat, const vec4_x ind)
{
	val_x temp[4*5];
	int i,j,k;
	int lInd[4] = {0,1,2,3}; // row indices
	int pivotIndex;
	
	for(i=0; i<4; i++)
		for(j=0; j<4; j++)
			temp[i*5+j] = mat[i*4+j];
	for(j=0; j<4; j++)
		temp[j*5+4] = ind[j]; // last column
	
	// Gauss descent
	for(k=0; k<3; k++)
	{
		// row switch 
		pivotIndex = k;
		for(i=k+1; i<4; i++)
		{
			if(fabs_x(temp[lInd[i]*5+k]) > fabs_x(temp[lInd[pivotIndex]*5+k]))
			{
				pivotIndex = i;
			}
			
		}
		//if(EqualZero(temp[pivotIndex*4+k])) 
		if(fabs_x(temp[lInd[pivotIndex]*5+k]) == 0.0)
		{    
			//NSLog([self description]); // got some !?!
			return false;
		}
		else
		{
			int rr = lInd[pivotIndex];
			
			lInd[pivotIndex] = lInd[k];
			lInd[k] = rr;
		}
		
		// current column elimination
		for(i=k+1; i<4; i++)
		{
			temp[lInd[i]*5+k] = temp[lInd[i]*5+k] / temp[lInd[k]*5+k]; // lik calculation
			
			for(j=k+1; j<5; j++)
			{
				temp[lInd[i]*5+j] = temp[lInd[i]*5+j] - temp[lInd[i]*5+k]*temp[lInd[k]*5+j];
			}
		}
	}
	
	if(fabs_x(temp[lInd[3]*5+3]) == 0.0) 
	{    
		//NSLog([self description]); // did not get
		return false;
	}
	else
	{
		temp[lInd[3]*5+4] = temp[lInd[3]*5+4] / temp[lInd[3]*5+3];
		
		// big bug here ? check that it works for 4x4 matrices (was i=3; i>=...)
		for(i=2; i>=0; i--)
		{
			val_x sum = 0.0;
			
			if(fabs_x(temp[lInd[i]*5+i]) == 0.0) 
			{    
				//NSLog([self description]); // did not get
				return false;
			}
			for(j=i+1; j<4; j++)
			{
				sum += temp[lInd[i]*5+j] * temp[lInd[j]*5+4];
			}
			temp[lInd[i]*5+4] = (temp[lInd[i]*5+4] - sum) / temp[lInd[i]*5+i];
		}
	}
	
	res[0] = temp[lInd[0]*5+4];
	res[1] = temp[lInd[1]*5+4];
	res[2] = temp[lInd[2]*5+4];
	res[3] = temp[lInd[3]*5+4];
	
	return true;
}


#pragma mark mat4d
matn_x c3fx(mat4_set_identity)(mat4_x res)
{
	int i,j;
	
	for(i=0; i<4; i++)
		for(j=0; j<4; j++)
		{
			res[MIJ4(i,j)] = 0.0;
		}
	
	for(j=0; j<4; j++)
	{
		res[MIJ4(j,j)] = 1.0;
	}
	
	return res;
}
matn_x c3fx(mat4_set_translation)(mat4_x res, val_x x, val_x y, val_x z)
{
	c3fx(mat4_set_identity)(res);
	res[MIJ4(0,3)] = x;
	res[MIJ4(1,3)] = y;
	res[MIJ4(2,3)] = z;
	
	return res;
}

matn_x c3fx(mat4_set_scale)(mat4_x res, val_x xs, val_x ys, val_x zs)
{
	c3fx(mat4_set_identity)(res);
	res[MIJ4(0,0)] = xs;
	res[MIJ4(1,1)] = ys;
	res[MIJ4(2,2)] = zs;
	
	return res;
}
matn_x c3fx(mat4_set_rotation)(mat4_x res, val_x x, val_x y, val_x z)
{
	// rotx is applied first, then y, then z.
	val_x cx = cos_x(x), cy = cos_x(y), cz = cos_x(z);
	val_x sx = sin_x(x), sy = sin_x(y), sz = sin_x(z);
	
	res[MIJ4(0,0)] = cy*cz;
	res[MIJ4(0,1)] = sx*sy*cz-cx*sz;
	res[MIJ4(0,2)] = cx*sy*cz+sx*sz;
	res[MIJ4(0,3)] = 0.0;
	
	res[MIJ4(1,0)] = cy*sz;
	res[MIJ4(1,1)] = sx*sy*sz+cx*cz;
	res[MIJ4(1,2)] = cx*sy*sz-sx*cz;
	res[MIJ4(1,3)] = 0.0;
	
	res[MIJ4(2,0)] = -sy;
	res[MIJ4(2,1)] = sx*cy;
	res[MIJ4(2,2)] = cx*cy;
	res[MIJ4(2,3)] = 0.0;
	
	res[MIJ4(3,0)] = 0.0;
	res[MIJ4(3,1)] = 0.0;
	res[MIJ4(3,2)] = 0.0;
	res[MIJ4(3,3)] = 1.0;
	
	return res;
}
matn_x c3fx(mat4_set_rotation_quaternion)(mat4_x res, const quat4_x q)
{
	val_x w = q[0], x = q[1], y = q[2], z = q[3];
	val_x x2 = C3SQR(q[1]), y2 = C3SQR(q[2]), z2 = C3SQR(q[3]);
	
	res[MIJ4(0,3)] = 0.0;
	res[MIJ4(1,3)] = 0.0;
	res[MIJ4(2,3)] = 0.0;
	res[MIJ4(3,3)] = 1.0;
	res[MIJ4(3,0)] = 0.0;
	res[MIJ4(3,1)] = 0.0;
	res[MIJ4(3,2)] = 0.0;
	
	res[MIJ4(0,0)] = 1.0-2.0*(y2+z2);
	res[MIJ4(1,1)] = 1.0-2.0*(x2+z2);
	res[MIJ4(2,2)] = 1.0-2.0*(x2+y2);
	
	res[MIJ4(0,1)] = 2.0*(x*y - w*z);
	res[MIJ4(0,2)] = 2.0*(w*y + x*z);
	res[MIJ4(1,0)] = 2.0*(x*y + w*z);
	res[MIJ4(1,2)] = 2.0*(y*z - w*x);
	res[MIJ4(2,0)] = 2.0*(x*z - w*y);
	res[MIJ4(2,1)] = 2.0*(y*z + w*x);	
	
	return res;
}
matn_x c3fx(mat4_set_rotation_rodriguez)(mat4_x res, const vec3_x axis_ampl)
{
#define RODRIG_EPSILON 0.00001
	vec3_x axis;
	val_x theta = c3fx(vec3_norm)(axis_ampl);
	
	c3fx(mat4_set_identity)(res);
	
	if(theta < RODRIG_EPSILON) 
	{
		return res;
	}
	else
	{
		c3fx(vec3_scale)(axis, axis_ampl, 1./theta);
		
		return c3fx(mat4_set_rotation_axis_angle)(res, axis, theta);
	}
}
matn_x c3fx(mat4_set_transform)(mat4_x res, const vec3_x t, const vec3_x r, const vec3_x s)
{
	mat4_x tm, rm, sm, tmp;
	
	c3fx(mat4_set_translation)(tm, t[0], t[1], t[2]);
	c3fx(mat4_set_rotation)(rm, r[0], r[1], r[2]);
	c3fx(mat4_set_scale)(sm, s[0], s[1], s[2]);
	
	c3fx(mat4_multiply)(tmp, r, s);
	return c3fx(mat4_multiply)(res, t, tmp);
}
void c3fx(mat4_get_transform)(const mat4_x m, vec3_x t, vec3_x r, vec3_x s)
{
	t[0] = m[MIJ4(0,3)];
	t[1] = m[MIJ4(1,3)];
	t[2] = m[MIJ4(2,3)];
	
	s[0] = sqrt_x(C3SQR(m[MIJ4(0,0)]) + C3SQR(m[MIJ4(1,0)]) + C3SQR(m[MIJ4(2,0)]));
	s[1] = sqrt_x(C3SQR(m[MIJ4(0,1)]) + C3SQR(m[MIJ4(1,1)]) + C3SQR(m[MIJ4(2,1)]));
	s[2] = sqrt_x(C3SQR(m[MIJ4(0,2)]) + C3SQR(m[MIJ4(1,2)]) + C3SQR(m[MIJ4(2,2)]));
	
	r[0] = 0.0;
	r[1] = 0.0;
	r[2] = 0.0;
	
	if(fabs_x(s[0]) > 0.0 && fabs_x(s[1]) > 0.0 && fabs_x(s[2]) > 0.0)
	{
		double r20 = m[MIJ4(2,0)] / s[0];
		
		if(r20 > 0.99999) // gimbal lock, X and Z axis are the same because of Y rotation
		{
			r[1] = -M_PI*0.5;
			r[0] = 0.0;
			r[2] = atan2_x(m[MIJ4(0,1)], m[MIJ4(1,1)]);
		}
		else if(r20 < -0.99999)
		{
			r[1] = M_PI*0.5;
			r[0] = 0.0;
			r[2] = atan2_x(m[MIJ4(0,1)], m[MIJ4(1,1)]);
		}
		else
		{
			r[1] = -asin_x(m[MIJ4(2,0)] / s[0]);
			r[0] = atan2_x(m[MIJ4(2,1)] / s[1], m[MIJ4(2,2)] / s[2]);
			r[2] = atan2_x(m[MIJ4(1,0)], m[MIJ4(0,0)]);
		}
	}
}
c3extern matn_x c3fx(mat4_set_look_at)(mat4_x res, const vec3_x eye, const vec3_x target, const vec3_x up)
{
	vec3_x tx, ty, tz;
	
	c3fx(vec3_subtract)(tz, target, eye);
	c3fx(vec3_normalize)(tz, tz);
	c3fx(vec3_cross_product)(tx, tz, up);
	c3fx(vec3_normalize)(tx, tx);
	c3fx(vec3_cross_product)(ty, tx, tz);
	c3fx(vec3_normalize)(ty, ty);
	c3fx(vec3_negate)(tz, tz);
	
	c3fx(mat4_set_identity)(res);
	
	// set these axes as column for the matrix
	res[0] = tx[0];
	res[4] = tx[1];
	res[8] = tx[2];
	
	res[1] = ty[0];
	res[5] = ty[1];
	res[9] = ty[2];
	
	res[2] = tz[0];
	res[6] = tz[1];
	res[10] = tz[2];
	
	res[3] = eye[0];
	res[7] = eye[1];
	res[11] = eye[2];
	
	return res;
}
matn_x c3fx(mat4_set_rotation_axis_angle)(mat4_x res, const normal_x axis, val_x angle)
{
	val_x xy, xz, yz, xs, ys, zs, s, c, t;
	
	s = sin_x(angle);
	c = cos_x(angle);
	t = 1.0 - c;
	
	xy = axis[0] * axis[1];
	xz = axis[0] * axis[2];
	yz = axis[1] * axis[2];
	xs = - s * axis[0];
	ys = - s * axis[1];
	zs = - s * axis[2];
	
	res[MIJ4(0,0)] = t * axis[0] * axis[0] + c;
	res[MIJ4(0,1)] = t * xy + zs;
	res[MIJ4(0,2)] = t * xz - ys;
	
	res[MIJ4(1,0)] = t * xy - zs;
	res[MIJ4(1,1)] = t * axis[1] * axis[1] + c ;
	res[MIJ4(1,2)] = t * yz + xs;
	
	res[MIJ4(2,0)] = t * xz + ys;
	res[MIJ4(2,1)] = t * yz - xs;
	res[MIJ4(2,2)] = t * axis[2] * axis[2] + c ;
	
	res[MIJ4(3,0)] = 0.0; res[MIJ4(3,1)] = 0.0; res[MIJ4(3,2)] = 0.0; 
	res[MIJ4(0,3)] = 0.0; res[MIJ4(1,3)] = 0.0; res[MIJ4(2,3)] = 0.0;
	res[MIJ4(3,3)] = 1.0;
	
	return res;
}

matn_x c3fx(mat4_set_ortho)(mat4_x res, val_x left, val_x right, val_x bottom, val_x top, val_x near, val_x far)
{
	c3fx(mat4_set_identity)(res);
	
	res[MIJ4(0,0)] = 2.0 / (right-left);
	res[MIJ4(1,1)] = 2.0 / (top-bottom);
	res[MIJ4(2,2)] = -2.0 / (far-near);
	res[MIJ4(0,3)] = - (right+left) / (right-left);
	res[MIJ4(1,3)] = - (top+bottom) / (top-bottom);
	res[MIJ4(2,3)] = - (far+near)   / (far-near);
	
	return res;
}
matn_x c3fx(mat4_set_perspective)(mat4_x res, val_x left, val_x right, val_x bottom, val_x top, val_x near, val_x far)
{
	c3fx(mat4_set_identity)(res);
	
	res[MIJ4(0,0)] = 2.0*near / (right-left);
	res[MIJ4(1,1)] = 2.0*near / (top-bottom);
	res[MIJ4(2,2)] = - (far+near) / (far-near);
	res[MIJ4(0,2)] = (right+left) / (right-left);
	res[MIJ4(1,2)] = (top+bottom) / (top-bottom);
	res[MIJ4(2,3)] = - 2.0*(far*near)   / (far-near);
	res[MIJ4(3,2)] = -1.0;
	res[MIJ4(3,3)] = 0.0;
	
	return res;
}
matn_x c3fx(mat4_set_perspective_fovx)(mat4_x res, val_x fovx, val_x aspect, val_x near, val_x far)
{
	val_x r = near * tan_x(fovx*0.5);
	
	return c3fx(mat4_set_perspective)(res, -r, r, -r/aspect, r/aspect, near, far);
}

matn_x c3fx(mat4_set_perspective_fovy)(mat4_x res, val_x fovy, val_x aspect, val_x near, val_x far)
{
	val_x l = near * tan_x(fovy*0.5);
	
	return c3fx(mat4_set_perspective)(res, -l*aspect, l*aspect, -l, l, near, far);
}
matn_x c3fx(mat4_set_row)(mat4_x res, int index, const vec4_x r)
{
	res[MIJ4(index, 0)] = r[0];
	res[MIJ4(index, 1)] = r[1];
	res[MIJ4(index, 2)] = r[2];
	res[MIJ4(index, 3)] = r[3];
	
	return res;
}
matn_x c3fx(mat4_set_column)(mat4_x res, int index, const vec4_x c)
{
	res[MIJ4(0, index)] = c[0];
	res[MIJ4(1, index)] = c[1];
	res[MIJ4(2, index)] = c[2];
	res[MIJ4(3, index)] = c[3];
	
	return res;
}
vecn_x c3fx(mat4_get_row)(vec4_x res, const mat4_x m, int index)
{
	int i;
	
	for(i=0; i<4; i++)
	{
		res[i] = m[MIJ4(index,i)];
	}
	
	return res;
}
vecn_x c3fx(mat4_get_column)(vec4_x res, const mat4_x m, int index)
{
	int i;
	
	for(i=0; i<4; i++)
	{
		res[i] = m[MIJ4(i,index)];
	}
	
	return res;
}
bool c3fx(mat4_compare)(const mat4_x m1, const mat4_x m2, val_x tol)
{
	int i;
	
	for(i=0;i<16; i++) if (fabs_x(m1[i]-m2[i]) > tol) return false;
	
	return true;
}
bool c3fx(mat4_equal)(const mat4_x m1, const mat4_x m2)
{
	return c3fx(mat4_compare)(m1, m2, 0.0);
}
bool c3fx(mat4_is_identity)(const mat4_x m1)
{
	mat4_x ident;
	
	return c3fx(mat4_equal)(m1, c3fx(mat4_set_identity)(ident));
}
matn_x c3fx(mat4_transpose)(mat4_x res, const mat4_x src)
{
	int i,j;
	
	for(i=0;i<4;i++)
	{
		for(j=0;j<4;j++)
		{
			res[MIJ4(i,j)] = src[MIJ4(j,i)];
		}
	}
	
	return res;
}
bool c3fx(mat4_invert)(mat4_x res, const mat4_x src)
{
	assert(res != src);
	
	vec4_x vec = {1.0, 0.0, 0.0, 0.0};
	vec4_x sol;
	
	if(!c3fx(vec4_solve_system)(sol, src, vec)) return false;
	c3fx(mat4_set_column)(res, 0, sol);
	
	vec[0] = 0.0; vec[1] = 1.0;
	if(!c3fx(vec4_solve_system)(sol, src, vec)) return false;
	c3fx(mat4_set_column)(res, 1, sol);
	
	vec[1] = 0.0; vec[2] = 1.0;
	if(!c3fx(vec4_solve_system)(sol, src, vec)) return false;
	c3fx(mat4_set_column)(res, 2, sol);
	
	vec[2] = 0.0; vec[3] = 1.0;
	if(!c3fx(vec4_solve_system)(sol, src, vec)) return false;
	c3fx(mat4_set_column)(res, 3, sol);
	
	return true;
	
}
bool c3fx(mat4_invert_transpose)(mat4_x res, const mat4_x src)
{
	assert(res != src);
	
	mat4_x tmp;
	
	return c3fx(mat4_invert)(res, c3fx(mat4_transpose)(tmp, src));
}
void    c3fx(mat4_print)(const char *prefix, const mat4_x _matrix)
{
	printf("%s\n[%f %f %f %f]\n[%f %f %f %f]\n[%f %f %f %f]\n[%f %f %f %f]", prefix, 
	       _matrix[0], _matrix[1], _matrix[2], _matrix[3], _matrix[4], _matrix[5], _matrix[6], _matrix[7], 
	       _matrix[8], _matrix[9], _matrix[10], _matrix[11], _matrix[12], _matrix[13], _matrix[14], _matrix[15]);
}


#pragma mark ray

bool c3fx(ray_intersect_plane)(vec3_x res, const ray_x ray, const plane_x plane)
{
	val_x den;
	vec3_x dirCoords;
	
	c3fx(ray_get_direction)(dirCoords, ray);
	
	den = plane[0]*dirCoords[0] + plane[1]*dirCoords[1] + plane[2]*dirCoords[2];
	
	if(den == 0.0) return false;
	else
	{
		val_x num, alpha;
		
		num = plane[0]*ray[0] + plane[1]*ray[1] + plane[2]*ray[2] + plane[3];
		alpha = -num/den;
		res[0] = ray[0] + alpha*dirCoords[0]; 
		res[1] = ray[1] + alpha*dirCoords[1]; 
		res[2] = ray[2] + alpha*dirCoords[2]; 
		
		return true;
	}
}
bool c3fx(ray_intersect_plane_ex)(vec3_x res, val_x *t, const ray_x ray, const plane_x plane)
{
	val_x den;
	vec3_x dirCoords;
	
	c3fx(ray_get_direction)(dirCoords, ray);
	
	den = plane[0]*dirCoords[0] + plane[1]*dirCoords[1] + plane[2]*dirCoords[2];
	
	if(den == 0.0) return false;
	else
	{
		val_x num, alpha;
		
		num = plane[0]*ray[0] + plane[1]*ray[1] + plane[2]*ray[2] + plane[3];
		alpha = -num/den;
		res[0] = ray[0] + alpha*dirCoords[0]; 
		res[1] = ray[1] + alpha*dirCoords[1]; 
		res[2] = ray[2] + alpha*dirCoords[2];  
		*t = alpha;
		
		return true;
	}
}
static int c3fx(line_intersect_segment)(const val_x p[2], const val_x a[2], const val_x b[2])
{
	val_x u1, u2, v1, v2, vdiff, t;
	
	u1 = a[0] - p[0]; v1 = a[1] - p[1];
	u2 = b[0] - p[0]; v2 = b[1] - p[1];
	
	if((v1 > 0.0 && v2 > 0.0) || (v1 < 0.0 && v2 < 0.0) || (u1 < 0.0 && u2 < 0.0))
		return 0;
	
	if(u1 > 0.0 && u2 > 0.0) return 1;
	
	vdiff = v2 - v1;
	if(fabs_x(vdiff) < 0.000001) 
	{
		if(fabs_x(v1) > 0.000001 || u1 > 0.0 || u2 > 0.0) return 0;
		return 1;
	}
	t = -v1 / vdiff;
	return ((u1 + t * (u2-u1) ) > 0.0)?1:0;
}
bool c3fx(ray_intersect_triangle)(vec3_x res, val_x *pt, val_x mint, val_x maxt, const vec3_x at, const vec3_x bt, const vec3_x ct, 
							 const normal_x normal, const ray_x ray)
{
	plane_x plane;
	vec3_x currentIntersection;
	val_x t;
	
	c3fx(plane_set_normal_point)(plane, normal, at);
	
	if(c3fx(ray_intersect_plane_ex)(currentIntersection, &t, ray, plane))
	{
		*pt = t;
		c3fx(vec3_copy)(res, currentIntersection);
		if(t >= mint && t < maxt)
		{
			int axis;
			int intCount;
			val_x a[2], b[2], c[2], p[2], minp[2], maxp[3], *n = (val_x*)normal;
			
			
			if(fabs_x(n[0]) >= fabs_x(n[1]) && fabs_x(n[0]) >= fabs_x(n[2])) axis = 0;
			else if(fabs_x(n[1]) >= fabs_x(n[0]) && fabs_x(n[1]) >= fabs_x(n[2])) axis = 1;
			else axis = 2;
			
			// we choose the best axis
			a[0] = at[(axis+1)%3];
			a[1] = at[(axis+2)%3];
			b[0] = bt[(axis+1)%3];
			b[1] = bt[(axis+2)%3];
			c[0] = ct[(axis+1)%3];
			c[1] = ct[(axis+2)%3];
			p[0] = currentIntersection[(axis+1)%3];
			p[1] = currentIntersection[(axis+2)%3];
			
			// computing triangle bounds
			minp[0] = C3MIN(a[0], C3MIN(b[0], c[0]));
			maxp[0] = C3MAX(a[0], C3MAX(b[0], c[0]));
			minp[1] = C3MIN(a[1], C3MIN(b[1], c[1]));
			maxp[1] = C3MAX(a[1], C3MAX(b[1], c[1]));
			
			// triangle bounds check
			if(p[0] < minp[0] || p[1] < minp[1] || p[0] > maxp[0] || p[1] > maxp[1])
				return false;
			
			intCount = 0;
			intCount += c3fx(line_intersect_segment)(p, a, b);
			intCount += c3fx(line_intersect_segment)(p, b, c);
			intCount += c3fx(line_intersect_segment)(p, c, a);
			
			if(intCount & 1) // odd winding rule (but not relevent because of triangles !)
			{
				return true;
			}
		}
	}
	return false;
}
bool c3fx(ray_intersect_triangle_major_bounds)(vec3_x res, val_x *pt, val_x mint, val_x maxt, const vec3_x at, const vec3_x bt, const vec3_x ct, 
										  const normal_x normal, unsigned char major, const val_x bounds[4], const ray_x ray)
{
	plane_x plane;
	vec3_x currentIntersection;
	val_x t;
	
	c3fx(plane_set_normal_point)(plane, normal, at);
	
	if(c3fx(ray_intersect_plane_ex)(currentIntersection, &t, ray, plane))
	{
		*pt = t;
		c3fx(vec3_copy)(res, currentIntersection);
		if(t >= mint && t < maxt)
		{
			int axis = major;
			int intCount;
			val_x a[2], b[2], c[2], p[2];
			
			// we choose the best axis
			a[0] = at[(axis+1)%3];
			a[1] = at[(axis+2)%3];
			b[0] = bt[(axis+1)%3];
			b[1] = bt[(axis+2)%3];
			c[0] = ct[(axis+1)%3];
			c[1] = ct[(axis+2)%3];
			p[0] = currentIntersection[(axis+1)%3];
			p[1] = currentIntersection[(axis+2)%3];
			
			
			// triangle bounds check
			if(p[0] < bounds[0] || p[1] < bounds[1] || p[0] > bounds[2] || p[1] > bounds[3])
				return false;
			
			intCount = 0;
			intCount += c3fx(line_intersect_segment)(p, a, b);
			intCount += c3fx(line_intersect_segment)(p, b, c);
			intCount += c3fx(line_intersect_segment)(p, c, a);
			
			if(intCount & 1) // odd winding rule (but not relevent because of triangles !)
			{
				return true;
			}
		}
	}
	return false;
}
bool c3fx(ray_intersect_box)(const ray_x ray, const box_x box)
{
	val_x near, far;
	
	return c3fx(ray_intersect_box_ex)(&near, &far, ray, box);
}
bool c3fx(ray_intersect_box_ex)(val_x *near, val_x *far, const ray_x ray, const box_x box)
{
	int i;
	val_x t1, t2, di;
	
	*near = -MAXFLOAT; *far = MAXFLOAT;
	
	for(i=0; i<3; i++)
	{
		if(ray[3+i] != 0.0)
		{
			di = 1.0 / ray[3+i];
			t1 = (box[i] - ray[i]) * di;
			t2 = (box[3+i] - ray[i]) * di;
			
			if(t1 > t2) 
			{
				di = t1; t1 = t2; t2 = di;
			}
			
			if(t1 > *near) *near = t1;
			if(t2 < *far) *far = t2;
			if(*near > *far || (*far < 0.0)) return false;
		}
		else 
		{
			if(!(ray[i] >= box[i] && ray[i] <= box[3+i])) return false;
			else continue;
		}
	}
	/*
	 *near = tmin;
	 *far = tmax;
	 */
	return true;
}
vecn_x c3fx(ray_project_point)(vec3_x res, const ray_x ray, const vec3_x v)
{
	plane_x plane;
	
	c3fx(plane_set_normal_point)(plane, C3_RAY_DIRECTION(ray), v);
	c3fx(ray_intersect_plane)(res, ray, plane); // should always return true
	
	return res;
}
val_x c3fx(ray_distance_with_point)(const ray_x ray, const vec3_x p)
{
	return sqrt(c3fx(ray_square_distance_with_point)(ray, p));
}
val_x c3fx(ray_square_distance_with_point)(const ray_x ray, const vec3_x p)
{
	vec3_x v, cp;
	
	c3fx(vec3_subtract)(v, p, C3_RAY_ORIGIN(ray));
	c3fx(vec3_cross_product)(cp, C3_RAY_DIRECTION(ray), v);
	return c3fx(vec3_square_norm)(cp)/c3fx(vec3_square_norm)(C3_RAY_DIRECTION(ray));
}
vecn_x c3fx(ray_symmetric_point)(vec3_x res, const ray_x ray, const vec3_x src)
{
	vec3_x proj, diff;
	
	c3fx(ray_project_point)(proj, ray, src);
	c3fx(vec3_subtract)(diff, proj, src);
	
	return c3fx(vec3_add)(res, proj, diff);
}
rayn_x c3fx(ray_transform)(ray_x res, const mat4_x mat, const ray_x src)
{
	vec3_x to, te, o, e;
	
	c3fx(vec3_add)(e, C3_RAY_ORIGIN(src), C3_RAY_DIRECTION(src));
	c3fx(vec3_copy)(o, C3_RAY_ORIGIN(src));
	c3fx(vec3_transform)(to, mat, o);
	c3fx(vec3_transform)(te, mat, e);
	c3fx(vec3_subtract)(e, te, to);
	
	return c3fx(ray_set)(res, to, e);
}

/* returns in res the points on b closest to a */
bool c3fx(ray_point_nearest_to_ray)(vec3_x res, const ray_x b, const ray_x a)
{
	plane_x plane, tmp;
	
	c3fx(vec3_cross_product)(tmp, C3_RAY_DIRECTION(b), C3_RAY_DIRECTION(a));
	c3fx(vec3_cross_product)(plane, tmp, C3_RAY_DIRECTION(a));
	
	plane[3] = -c3fx(vec3_dot_product)(C3_RAY_ORIGIN(a), plane);
	return c3fx(ray_intersect_plane)(res, b, plane);
}

/* returns the oriented angle between ray/plane intersection and p, relatively to center */
val_x c3fx(ray_angle_with_plane_point_center)(const ray_x ray, const plane_x plane, const vec3_x center, const vec3_x p)
{
	vec3_x rref, rint;
	val_x a;
	
	c3fx(vec3_subtract)(rref, p, center);
	c3fx(ray_intersect_plane)(rint, ray, plane);
	c3fx(vec3_subtract)(rint, rint, center);
	
	a = c3fx(vec3_oriented_angle)(rint, rref, plane);
	if(a >= M_PI) a -= 2.0*M_PI;
	
	return a;
}
vecn_x c3fx(triangle_get_normal)(vec3_x res, const vec3_x a, const vec3_x b, const vec3_x c)
{
	vec3_x ab, bc;
	
	c3fx(vec3_subtract)(ab, b, a);
	c3fx(vec3_subtract)(bc, c, b);
	c3fx(vec3_cross_product)(res, ab, bc);
	c3fx(vec3_normalize)(res, res);
	
	return res;
}

#pragma mark box
boxn_x c3fx(box_copy)(box_x res, const box_x src)
{
	res[0] = src[0];
	res[1] = src[1];
	res[2] = src[2];
	res[3] = src[3];
	res[4] = src[4];
	res[5] = src[5];
	
	return res;
}
boxn_x c3fx(box_set_null)(box_x res)
{
	res[0] = 0.0;
	res[1] = 0.0;
	res[2] = 0.0;
	res[3] = 0.0;
	res[4] = 0.0;
	res[5] = 0.0;
	
	return res;
}
boxn_x c3fx(box_set)(box_x res, const vec3_x low, const vec3_x high)
{
	c3fx(vec3_copy)(&res[0], low);
	c3fx(vec3_copy)(&res[3], high);
	
	return res;
}
boxn_x c3fx(box_set_center_size)(box_x res, const vec3_x c, val_x xw, val_x yw, val_x zw)
{
	xw *= 0.5; yw *= 0.5; zw *= 0.5;
	
	res[0] = c[0] - xw;
	res[1] = c[1] - yw;
	res[2] = c[2] - zw;
	res[3] = c[0] + xw;
	res[4] = c[1] + yw;
	res[5] = c[2] + zw;
	
	return res;
}
boxn_x c3fx(box_set_from_triangle)(box_x res, const vec3_x a, const vec3_x b, const vec3_x c)
{
	res[0] = C3MIN(a[0], C3MIN(b[0], c[0]));
	res[1] = C3MIN(a[1], C3MIN(b[1], c[1]));
	res[2] = C3MIN(a[2], C3MIN(b[2], c[2]));
	
	res[3] = C3MAX(a[0], C3MAX(b[0], c[0]));
	res[4] = C3MAX(a[1], C3MAX(b[1], c[1]));
	res[5] = C3MAX(a[2], C3MAX(b[2], c[2]));
	
	return res;
}
boxn_x c3fx(box_union)(box_x res, const box_x a, const box_x b)
{
	if(c3fx(box_is_empty)(b)) return c3fx(box_copy)(res, a);
	if(c3fx(box_is_empty)(a)) return c3fx(box_copy)(res, b);
	
	res[0] = (a[0] < b[0] ? a[0] : b[0]);
	res[1] = (a[1] < b[1] ? a[1] : b[1]);
	res[2] = (a[2] < b[2] ? a[2] : b[2]);
	
	res[0+3] = (a[0+3] > b[0+3] ? a[0+3] : b[0+3]);
	res[1+3] = (a[1+3] > b[1+3] ? a[1+3] : b[1+3]);
	res[2+3] = (a[2+3] > b[2+3] ? a[2+3] : b[2+3]);
	
	return res;
}
boxn_x c3fx(box_intersection)(box_x res, const box_x a, const box_x b)
{
	if(c3fx(box_is_empty)(b)) return c3fx(box_copy)(res, b);
	if(c3fx(box_is_empty)(a)) return c3fx(box_copy)(res, a);
	
	res[0] = C3MAX(a[0], b[0]);
	res[1] = C3MAX(a[1], b[2]);
	res[2] = C3MAX(a[2], b[2]);
	
	res[3] = C3MIN(a[3], b[3]);
	res[4] = C3MIN(a[4], b[4]);
	res[5] = C3MIN(a[5], b[5]); 
	
	if(res[0] >= res[3] || res[1] >= res[4] || res[2] >= res[5])
	{
		c3fx(box_set_null)(res);
	}
	
	return res;
}
boxn_x c3fx(box_add_point)(box_x res, const box_x src, vec3_x p)
{
	if(c3fx(box_is_empty)(src)) return c3fx(box_copy)(res, src);
	
	res[0] = (src[0] < p[0] ? src[0] : p[0]);
	res[1] = (src[1] < p[1] ? src[1] : p[1]);
	res[2] = (src[2] < p[2] ? src[2] : p[2]);
	
	res[0+3] = (src[0+3] > p[0] ? src[0+3] : p[0]);
	res[1+3] = (src[1+3] > p[1] ? src[1+3] : p[1]);
	res[2+3] = (src[2+3] > p[2] ? src[2+3] : p[2]);
	
	return res;
}
bool c3fx(box_is_empty)(const box_x box)
{
	return box[0] >= box[3] || box[1] >= box[4] || box[2] >= box[5];
}
bool c3fx(box_instersect_sphere)(const box_x box, const vec3_x center, val_x radius)
{
	val_x dmax = 0.0, dmin = 0.0, a, b, r2 = radius*radius;
	int i;
	
	for(i=0; i<3; i++)
	{
		a = center[i]-box[i]; a=a*a;
		b = center[i]-box[3+i]; b=b*b;
		dmax += C3MAX(a,b);
		
		if(center[i] < box[i]) dmin += a;
		else if(center[i] > box[3+i]) dmin += b;
	}
	
	return (r2 >= dmin && r2 <= dmax);
}
bool c3fx(box_contain_point)(const box_x box, const vec3_x p)
{
	if(p[0] < box[0] || p[0] > box[3+0]
	   || p[1] < box[1] || p[1] > box[3+1]
	   || p[2] < box[2] || p[2] > box[3+2]) return false;
	
	return true;
}
boxn_x c3fx(box_transform)(box_x res, const mat4_x mat, const box_x src)
{
	vec3_x vo, vi;
	int i=0;
	
	res[0] = res[1] = res[2] = MAXFLOAT;
	res[3] = res[4] = res[5] = -MAXFLOAT;
	
	for(i=0; i<8; i++)
	{
		vo[0] = (i&1)?src[3]:src[0];
		vo[1] = (i&2)?src[4]:src[1];
		vo[2] = (i&4)?src[5]:src[2];
		
		c3fx(vec3_transform)(vi, mat, vo);
		res[0] = C3MIN(vi[0], res[0]);
		res[1] = C3MIN(vi[1], res[1]);
		res[2] = C3MIN(vi[2], res[2]);
		
		res[3] = C3MAX(vi[0], res[3]);
		res[4] = C3MAX(vi[1], res[4]);
		res[5] = C3MAX(vi[2], res[5]);
	}
	
	return res;
}

#include "C3TriangleBoxOverlap-generic.h"
//extern int _TDTriangleOverlapsBox(double boxcenter[3],double boxhalfsize[3],double triverts[3][3]);

bool c3fx(box_instersect_triangle)(const box_x box, const vec3_x a, const vec3_x b, const vec3_x c)
{
	val_x cb[3], sb[3], vert[3][3];
	
	sb[0] = (box[3] - box[0]) * 0.5;
	sb[1] = (box[4] - box[1]) * 0.5;
	sb[2] = (box[5] - box[2]) * 0.5;
	
	cb[0] = (box[3] + box[0]) * 0.5;
	cb[1] = (box[4] + box[1]) * 0.5;
	cb[2] = (box[5] + box[2]) * 0.5; 
	
	vert[0][0] = a[0]; vert[0][1] = a[1]; vert[0][2] = a[2];
	vert[1][0] = b[0]; vert[1][1] = b[1]; vert[1][2] = b[2];
	vert[2][0] = c[0]; vert[2][1] = c[1]; vert[2][2] = c[2];
	
	return c3fx(_TDTriangleOverlapsBox)(cb, sb, vert);
}
#pragma mark frustum

frustumn_x c3fx(frustum_set_perspective)(frustum_x res, val_x hfovrad, rect_x norm_rect, val_x near, val_x far)
{
	val_x ta2 = tan(hfovrad*0.5);
	
	for(int i=0; i<4; i++)
	{
		vec2_x p;
		
		c3fx(rect_point_at_index)(p, norm_rect, i);
		
		res[i*3+0] = p[0] * near * ta2;
		res[i*3+1] = p[1] * near * ta2;	
		res[i*3+2] = -near;
		
		res[(i+4)*3+0] = p[0] * far * ta2;
		res[(i+4)*3+1] = p[1] * far * ta2;
		res[(i+4)*3+2] = -far;
	}
	
	return res;
}
frustumn_x c3fx(frustum_set_ortho)(frustum_x res, val_x fovwidth, rect_x norm_rect, val_x near, val_x far)
{
	val_x w2 = fovwidth * 0.5;
	
	for(int i=0; i<4; i++)
	{
		vec2_x p;
		
		c3fx(rect_point_at_index)(p, norm_rect, i);
		
		res[i*3+0] = p[0] * w2;
		res[i*3+1] = p[1] * w2;	
		res[i*3+2] = -near;
		
		res[(i+4)*3+0] = p[0] * w2;
		res[(i+4)*3+1] = p[1] * w2;
		res[(i+4)*3+2] = -far;
	}
	
	return res;
}
frustumn_x c3fx(frustum_transform)(frustum_x res, const mat4_x mat, const frustum_x v)
{
	for(int i=0; i<8; i++)
	{
		c3fx(vec3_transform)(&res[i*3], mat, &v[i*3]);
		
		// generic tests
		//vec4_copy_d(res, res);
		//float test[12];
		//vec4_copy_f(test, test);
		
	}
	
	return res;
}
bool c3fx(frustum_perspective_intersect_triangle)(const frustum_x frustum, const vec3_x at, const vec3_x bt, const vec3_x ct)
{
	vec3_x ab, bc, ca;	
	
	c3fx(vec3_subtract)(ab, at, bt);
	c3fx(vec3_subtract)(bc, bt, ct);
	c3fx(vec3_subtract)(ca, ct, at);
	
	return false;
}
bool c3fx(frustum_perspective_intersect_edge)(const frustum_x frustum, const vec3_x ae, const vec3_x be)
{
	vec3_x ab;	
	
	vec3_x e[6], a[11];
	
	// SAT algorithm: projection on face normals & edge cross products
	// get frustum edges
	c3fx(vec3_subtract)(e[0], frustum+3*1, frustum+3*0);
	c3fx(vec3_subtract)(e[1], frustum+3*3, frustum+3*0);
	c3fx(vec3_subtract)(e[2], frustum+3*4, frustum+3*0);
	c3fx(vec3_subtract)(e[3], frustum+3*5, frustum+3*1);
	c3fx(vec3_subtract)(e[4], frustum+3*6, frustum+3*2);
	c3fx(vec3_subtract)(e[5], frustum+3*7, frustum+3*3);
	
	// get frustum normals
	c3fx(vec3_cross_product)(a[6], e[0], e[1]);
	c3fx(vec3_cross_product)(a[7], e[0], e[2]);
	c3fx(vec3_cross_product)(a[8], e[1], e[3]); // because parallel
	c3fx(vec3_cross_product)(a[9], e[0], e[4]); // idem
	c3fx(vec3_cross_product)(a[10], e[1], e[5]);
	
	// edge cross products
	c3fx(vec3_subtract)(ab, ae, be);
	int i;
	for(i=0; i<6; i++)
		c3fx(vec3_cross_product)(a[i], e[i], ab);
	
	// now project along those 11 axis
	for(i=0; i<11; i++)
	{
		vec2_x frustum_bounds, edge_bounds;
		
		// normalize not needed, just for easier debugging
		//c3fx(vec3_normalize)(a[i], a[i]);
		c3fx(vec3_project_frustum)(frustum_bounds, a[i], frustum);
		c3fx(vec3_project_edge)(edge_bounds, a[i], ae, be);
		
		if(frustum_bounds[0] > edge_bounds[1] || edge_bounds[0] > frustum_bounds[1])
			return false;
	}
	
	return true;
}


