
// vec2
#pragma mark -
#pragma mark vec2
#pragma mark -
c3inline vecn_x c3fx(vec2_copy)(vec2_x res, const vec2_x src);
c3inline vecn_x c3fx(vec2_copy)(vec2_x res, const vec2_x src);
c3inline vecn_x c3fx(vec2_add)(vec2_x res, const vec2_x a, const vec2_x b);
c3inline vecn_x c3fx(vec2_subtract)(vec2_x res, const vec2_x a, const vec2_x b);
c3inline vecn_x c3fx(vec2_rotate)(vec2_x res, const vec2_x src, val_x angle);
c3inline vecn_x c3fx(vec2_lerp)(vec2_x res, const vec2_x a, const vec2_x b, val_x alpha);
c3inline val_x  c3fx(vec2_square_norm)(const vec2_x a);
c3inline val_x  c3fx(vec2_norm)(const vec2_x a);
c3inline val_x  c3fx(vec2_distance)(const vec2_x a, const vec2_x b);
c3inline vecn_x c3fx(vec2_multiply)(vec2_x res, const mat2_x mat, const vec2_x v);
c3inline vecn_x c3fx(vec2_transform)(vec2_x res, const mat3_x mat, const vec2_x v);

// mat2
#pragma mark -
#pragma mark mat2
#pragma mark -
c3inline matn_x c3fx(mat2_set_rotation)(mat2_x res, val_x a);
c3inline matn_x c3fx(mat2_set_scale)(mat2_x res, val_x sx, val_x sy);
c3inline matn_x c3fx(mat2_multiply)(mat2_x res, const mat2_x a, const mat2_x b);

// rect
#pragma mark -
#pragma mark rect
#pragma mark -
c3inline vecn_x c3fx(rect_point_at_index)(vec2_x res, const rect_x r, int index); 
c3inline vecn_x c3fx(rect_inset)(rect_x res, const rect_x r, val_x dx, val_x dy); 

// vec3
#pragma mark -
#pragma mark vec3
#pragma mark -
c3inline vecn_x c3fx(vec3_copy)(vec3_x res, const vec3_x src);
c3inline vecn_x c3fx(vec3_set)(vec3_x res, val_x x, val_x y, val_x z);
c3inline vecn_x c3fx(vec3_add)(vec3_x res, const vec3_x a, const vec3_x b);
c3inline vecn_x c3fx(vec3_subtract)(vec3_x res, const vec3_x a, const vec3_x b);
c3inline vecn_x c3fx(vec3_negate)(vec3_x res, const vec3_x src);
c3inline vecn_x c3fx(vec3_scale)(vec3_x res, const vec3_x src, val_x scalar);
c3inline bool   c3fx(vec3_compare)(const vec3_x a, const vec3_x b, val_x tol);
c3inline bool   c3fx(vec3_equal)(const vec3_x a, const vec3_x b);
c3inline bool   c3fx(vec3_is_null)(const vec3_x a);
c3inline val_x  c3fx(vec3_dot_product)(const vec3_x a, const vec3_x b);
c3inline vecn_x c3fx(vec3_cross_product)(vec3_x res, const vec3_x a, const vec3_x b);
c3inline val_x  c3fx(vec3_mixed_product)(const vec3_x a, const vec3_x b, const vec3_x c);
c3inline val_x  c3fx(vec3_square_norm)(const vec3_x a);
c3inline val_x  c3fx(vec3_norm)(const vec3_x a);
c3inline val_x  c3fx(vec3_distance)(const vec3_x a, const vec3_x b);
c3inline vecn_x c3fx(vec3_normalize)(vec3_x res, const vec3_x a);
c3inline vecn_x c3fx(vec3_homogenize)(vec3_x res, const vec3_x a);
c3inline vecn_x c3fx(vec3_add_scaled)(vec3_x res, const vec3_x p, const vec3_x q, val_x lambda);
c3inline vecn_x c3fx(vec3_lerp)(vec3_x res, const vec3_x a, const vec3_x b, val_x alpha);
c3inline vecn_x c3fx(vec3_multiply)(vec3_x res, const mat3_x mat, const vec3_x v);
c3inline vecn_x c3fx(vec3_transform)(vec3_x res, const mat4_x mat, const vec3_x v);
c3inline vecn_x c3fx(vec3_project_x)(vec3_x res, const vec3_x a);
c3inline vecn_x c3fx(vec3_project_y)(vec3_x res, const vec3_x a);
c3inline vecn_x c3fx(vec3_project_z)(vec3_x res, const vec3_x a);
c3inline vecn_x c3fx(vec3_project_xy)(vec3_x res, const vec3_x a);
c3inline vecn_x c3fx(vec3_project_yz)(vec3_x res, const vec3_x a);
c3inline vecn_x c3fx(vec3_project_zx)(vec3_x res, const vec3_x a);
c3inline void   c3fx(vec3_project_frustum)(vec2_x minmax, const vec3_x dir, const frustum_x frustum);
c3inline void   c3fx(vec3_project_edge)(vec2_x minmax, const vec3_x dir, const vec3_x at, const vec3_x bt);
c3inline void   c3fx(vec3_project_triangle)(vec2_x minmax, const vec3_x dir, const vec3_x at, const vec3_x bt, const vec3_x ct);
c3inline vecn_x c3fx(vec3_center_scale)(vec3_x res, const vec3_x a, const vec3_x center, val_x s);
c3inline vecn_x c3fx(vec3_reflect)(vec3_x res, const vec3_x a, const normal_x normal);
c3inline vecn_x c3fx(vec3_refract)(vec3_x res, const vec3_x a, const normal_x normal, val_x index);
c3extern vecn_x c3fx(vec3_rotate_centered_axis)(vec3_x res, const vec3_x src, val_x a, const normal_x axis, const vec3_x center);
c3extern val_x  c3fx(vec3_angle)(const vec3_x a, const vec3_x b);
c3extern val_x  c3fx(vec3_oriented_angle)(const vec3_x a, const vec3_x b, const vec3_x ref);

// normal
#pragma mark -
#pragma mark normal
#pragma mark -
c3inline vecn_x c3fx(normal_transform_it)(normal_x res, const mat4_x mat_it, const normal_x src);
c3inline bool   c3fx(normal_transform)(normal_x res, const mat4_x mat, const normal_x src);

// vec4
#pragma mark -
#pragma mark vec4
#pragma mark -
c3inline vecn_x c3fx(vec4_set)(vec4_x res, val_x x, val_x y, val_x z, val_x w);
c3inline vecn_x c3fx(vec4_copy)(vec4_x res, const vec4_x src);
c3inline bool   c3fx(vec4_compare)(const vec4_x a, const vec4_x b, val_x tol);
c3inline bool   c3fx(vec4_equal)(const vec4_x a, const vec4_x b);
c3inline vecn_x c3fx(vec4_lerp)(vec4_x res, const vec4_x a, const vec4_x b, val_x alpha);
c3inline vecn_x c3fx(vec4_homogenize)(vec4_x res, const vec4_x a);
c3inline vecn_x c3fx(vec4_multiply)(vec4_x res, const mat4_x mat, const vec4_x v);
c3inline vecn_x c3fx(vec4_multiply_delta)(vec4_x res, const mat4_x mat, const vec4_x v);
c3inline vecn_x c3fx(vec4_transpose_transform)(vec4_x res, const mat4_x m, const vec4_x v);
c3extern bool   c3fx(vec4_solve_system)(vec4_x res, const mat4_x mat, const vec4_x ind);

// quat4
#pragma mark -
#pragma mark quat4
#pragma mark -
c3inline quatn_x c3fx(quat4_set)(quat4_x res, val_x w, val_x x, val_x y, val_x z);
c3inline quatn_x c3fx(quat4_set_axis_angle)(quat4_x res, const normal_x axis, val_x angle);
c3inline quatn_x c3fx(quat4_set_from_rotation_matrix)(quat4_x res, const mat4_x rotmat);
c3inline quatn_x c3fx(quat4_copy)(quat4_x res, const quat4_x src);
c3inline quatn_x c3fx(quat4_product)(quat4_x res, const quat4_x q1, const quat4_x q2);
c3inline quatn_x c3fx(quat4_add)(quat4_x res, const quat4_x q1, const quat4_x q2);
c3inline quatn_x c3fx(quat4_subtract)(quat4_x res, const quat4_x q1, const quat4_x q2);
c3inline quatn_x c3fx(quat4_scale)(quat4_x res, const quat4_x q, val_x t);
c3inline val_x   c3fx(quat4_norm)(const quat4_x q);
c3inline quatn_x c3fx(quat4_normalize)(quat4_x res, const quat4_x q);
c3inline quatn_x c3fx(quat4_slerp)(quat4_x res, const quat4_x q1, const quat4_x q2, val_x t); 

// plane
#pragma mark -
#pragma mark plane
#pragma mark -
c3inline vecn_x c3fx(plane_copy)(plane_x res, const plane_x src);
c3inline vecn_x c3fx(plane_set_normal_point)(plane_x res, const normal_x normal, const vec3_x p);
c3inline vecn_x c3fx(plane_set_points)(plane_x res, const vec3_x a, const vec3_x b, const vec3_x c);
c3inline vecn_x c3fx(plane_transform_it)(plane_x res, const mat4_x mat_it, const plane_x src);
c3inline bool   c3fx(plane_transform)(plane_x res, const mat4_x mat, const plane_x src);


// mat4d
#pragma mark -
#pragma mark mat4
#pragma mark -
c3inline matn_x c3fx(mat4_copy)(mat4_x res, const mat4_x src);
c3inline matn_x c3fx(mat4_multiply)(mat4_x res, const mat4_x a, const mat4_x b);
c3extern matn_x c3fx(mat4_set_identity)(mat4_x res);
c3extern matn_x c3fx(mat4_set_identity)(mat4_x res);
c3extern matn_x c3fx(mat4_set_translation)(mat4_x res, val_x x, val_x y, val_x z);
c3extern matn_x c3fx(mat4_set_scale)(mat4_x res, val_x xs, val_x ys, val_x zs);
c3extern matn_x c3fx(mat4_set_rotation)(mat4_x res, val_x x, val_x y, val_x z);
c3extern matn_x c3fx(mat4_set_rotation_quaternion)(mat4_x res, const quat4_x q);
c3extern matn_x c3fx(mat4_set_rotation_axis_angle)(mat4_x res, const normal_x axis, val_x angle);
c3extern matn_x c3fx(mat4_set_rotation_rodriguez)(mat4_x res, const vec3_x axis_ampl);
c3extern matn_x c3fx(mat4_set_transform)(mat4_x res, const vec3_x t, const vec3_x r, const vec3_x s);
c3extern void   c3fx(mat4_get_transform)(const mat4_x m, vec3_x t, vec3_x r, vec3_x s);
c3extern matn_x c3fx(mat4_set_look_at)(mat4_x res, const vec3_x eye, const vec3_x target, const vec3_x up);
c3extern matn_x c3fx(mat4_set_ortho)(mat4_x res, val_x left, val_x right, val_x bottom, val_x top, val_x near, val_x far);

// near and far are expected to be positive, but the returned matrix is a -Z frustum
c3extern matn_x c3fx(mat4_set_perspective)(mat4_x res, val_x left, val_x right, val_x bottom, val_x top, val_x near, val_x far);
c3extern matn_x c3fx(mat4_set_perspective_fovx)(mat4_x res, val_x fovx, val_x aspect, val_x near, val_x far);
c3extern matn_x c3fx(mat4_set_perspective_fovy)(mat4_x res, val_x fovy, val_x aspect, val_x near, val_x far);

c3extern matn_x c3fx(mat4_set_row)(mat4_x res, int index, const vec4_x r);
c3extern matn_x c3fx(mat4_set_column)(mat4_x res, int index, const vec4_x c);
c3extern vecn_x c3fx(mat4_get_row)(vec4_x res, const mat4_x m, int index);
c3extern vecn_x c3fx(mat4_get_column)(vec4_x res, const mat4_x m, int index);
c3extern bool   c3fx(mat4_compare)(const mat4_x m1, const mat4_x m2, val_x tol);
c3extern bool   c3fx(mat4_equal)(const mat4_x m1, const mat4_x m2);
c3extern bool   c3fx(mat4_is_identity)(const mat4_x m);
c3extern matn_x c3fx(mat4_transpose)(mat4_x res, const mat4_x src);
c3extern bool   c3fx(mat4_invert)(mat4_x res, const mat4_x src);
c3extern bool   c3fx(mat4_invert_transpose)(mat4_x res, const mat4_x src);
c3extern void   c3fx(mat4_print)(const char *prefix, const mat4_x src);



// ray
#pragma mark -
#pragma mark ray
#pragma mark -
c3inline rayn_x c3fx(ray_set)(ray_x res, const vec3_x orig, const vec3_x dir);
c3inline rayn_x c3fx(ray_set_xaxis)(ray_x res);
c3inline rayn_x c3fx(ray_set_yaxis)(ray_x res);
c3inline rayn_x c3fx(ray_set_zaxis)(ray_x res);
c3inline rayn_x c3fx(ray_copy)(ray_x res, ray_x src);
c3inline vecn_x c3fx(ray_get_origin)(vec3_x res, const ray_x ray);
c3inline vecn_x c3fx(ray_get_direction)(vec3_x res, const ray_x ray);
c3extern bool   c3fx(ray_intersect_plane)(vec3_x res, const ray_x ray, const plane_x plane);
c3extern bool   c3fx(ray_intersect_plane_ex)(vec3_x res, val_x *t, const ray_x ray, const plane_x plane);
c3extern bool   c3fx(ray_intersect_triangle)(vec3_x res, val_x *pt, val_x mint, val_x maxt, const vec3_x at, const vec3_x bt, const vec3_x ct, 
										 const normal_x normal, const ray_x ray);
c3extern bool   c3fx(ray_intersect_triangle_major_bounds)(vec3_x res, val_x *pt, val_x mint, val_x maxt, const vec3_x at, 
													  const vec3_x bt, const vec3_x ct, 
													  const normal_x normal, unsigned char major, const val_x bounds[4], const ray_x ray);
c3extern bool   c3fx(ray_intersect_box)(const ray_x ray, const box_x box);
c3extern bool   c3fx(ray_intersect_box_ex)(val_x *near, val_x *far, const ray_x ray, const box_x box);
c3extern vecn_x c3fx(ray_project_point)(vec3_x res, const ray_x ray, const vec3_x v);
c3extern val_x  c3fx(ray_distance_with_point)(const ray_x ray, const vec3_x p);
c3extern val_x  c3fx(ray_square_distance_with_point)(const ray_x ray, const vec3_x p);
c3extern vecn_x c3fx(ray_symmetric_point)(vec3_x res, const ray_x ray, const vec3_x src);
c3extern rayn_x c3fx(ray_transform)(ray_x res, const mat4_x mat, const ray_x src);
c3extern bool   c3fx(ray_point_nearest_to_ray)(vec3_x res, const ray_x a, const ray_x b);// returns in res the points on a closest to b
c3extern val_x  c3fx(ray_angle_with_plane_point_center)(const ray_x ray, const plane_x plane, const vec3_x center, const vec3_x p);

// triangle
#pragma mark -
#pragma mark triangle
#pragma mark -
c3extern vecn_x c3fx(triangle_get_normal)(vec3_x res, const vec3_x a, const vec3_x b, const vec3_x c);

// box
#pragma mark -
#pragma mark box
#pragma mark -
c3extern boxn_x c3fx(box_copy)(box_x res, const box_x src);
c3extern boxn_x c3fx(box_set_null)(box_x res);
c3extern boxn_x c3fx(box_set)(box_x res, const vec3_x low, const vec3_x high);
c3extern boxn_x c3fx(box_set_center_size)(box_x res, const vec3_x c, val_x xw, val_x yw, val_x zw);
c3extern boxn_x c3fx(box_set_from_triangle)(box_x res, const vec3_x a, const vec3_x b, const vec3_x c);
c3extern boxn_x c3fx(box_union)(box_x res, const box_x a, const box_x b);
c3extern boxn_x c3fx(box_intersection)(box_x res, const box_x a, const box_x b);
c3extern boxn_x c3fx(box_add_point)(box_x res, const box_x src, vec3_x p);
c3extern bool   c3fx(box_is_empty)(const box_x box);
c3extern boxn_x c3fx(box_transform)(box_x res, const mat4_x mat, const box_x src);
c3extern bool   c3fx(box_contain_point)(const box_x box, const vec3_x p);
c3extern bool   c3fx(box_instersect_sphere)(const box_x box, const vec3_x center, val_x radius);
c3extern bool   c3fx(box_instersect_triangle)(const box_x box, const vec3_x a, const vec3_x b, const vec3_x c);

// frustum
#pragma mark -
#pragma mark frustum
#pragma mark -
c3extern frustumn_x c3fx(frustum_set_perspective)(frustum_x res, val_x hfovrad, rect_x norm_rect, val_x near, val_x far);
c3extern frustumn_x c3fx(frustum_set_ortho)(frustum_x res, val_x fovwidth, rect_x norm_rect, val_x near, val_x far);
c3extern frustumn_x c3fx(frustum_transform)(frustum_x res, const mat4_x mat, const frustum_x v);
c3extern bool    	c3fx(frustum_perspective_intersect_triangle)(const frustum_x frustum, const vec3_x at, const vec3_x bt, const vec3_x ct);
c3extern bool 		c3fx(frustum_perspective_intersect_edge)(const frustum_x frustum, const vec3_x ae, const vec3_x be);


