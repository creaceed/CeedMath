//
//  main.c
//  CeedMathDemo
//
//  Created by Raphael Sebbe on 04/12/11.
//  Copyright (c) 2011 Creaceed. All rights reserved.
//

#include <stdio.h>
#include <CeedMath/C3Math.h>

int main (int argc, const char * argv[])
{
	vec3_f a = {1,1,1}, b = {1,0,0}, c;
	val_f n;
	
	
	n = vec3_norm_f(a);
	printf("Norm of a: %f\n", n);
	
	vec3_subtract_f(c, a, b);
	n = vec3_norm_f(c);
	printf("Norm of c (c = a-b): %f\n", n);
	
	// chaining
	n = vec3_norm_f(vec3_subtract_f(c, a, b));
	
    return 0;
}

