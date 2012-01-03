About CeedMath
==============

CeedMath is a library for doing computer graphics math on Mac and iPhone, although it could be ported to other environment easily.

CeedMath is about computation with vectors and matrices, as well as boxes, rays, quaternions, and other CG types. It was designed to be easy to use and understand, and at the same time small and efficient.

CeedMath has simple data structures and operates easily with other formats. It is implemented in C with a form of template that allows a single implementation for both float and double types. That makes it possible to use it as is in Obj-C code. A polymorphic C++ version is also automatically derived.

Requirements
------------
CeedMath is mostly C based and should run on most compiler and tool versions. The provided Xcode project is made for iOS5 / Mac OS X 10.7.

Using CeedMath in your app
--------------------------
CeedMath project has been setup in a way that you can directly add it to your own application project. This enables easier updates of CeedMath with no need of separate build and file import (although you can still do it if that's what you want). 

Steps to include CeedMath:

* import the CeedMath project in your app project
* in your app target, add libCeedMath to the link phase (this creates an implicit build dependency)
* in your target build settings, add the following header search path (quotes matter):
	* `"$(CONFIGURATION_BUILD_DIR)/usr/local/include"`
	* `"$(DSTROOT)/usr/local/include"`
* you can then import CeedMath in your source file:
	* `#import <CeedMath/C3Math.h>`

A code example
--------------
Here's a C code example:

	vec3_f a = {1,1,1}, b = {1,0,0}, c;
	val_f n;
	
	
	n = vec3_norm_f(a);
	printf("Norm of a: %f\n", n);
	
	vec3_subtract_f(c, a, b);
	n = vec3_norm_f(c);
	printf("Norm of c (c = a-b): %f\n", n);

Note that in C++, you can even omit the `_f` or `_d` suffix for functions, because of the polymorphic implementation.


Why that name?
--------------
Well, it's how we name our reusable libs, both the internal and public ones here at Creaceed. English suffix "ceed" (as found in "succeed", "exceed", etc.) comes from Latin word "cedere", which means: to yield.

What's next?
------------
CeedMath is made to evolve. Feel free to add new methods, new types, or any other evolution. It is far from complete, but it provides a good starting point in the case you need CG math in your project.

Things that sound like possible in a not too distant future:

* In its current state, CeedMath has no support for SIMD (NEON, SSE…). This is also a possible evolution path for the library.
* Bridge to Apple's CG and NS types (CGPoint, NSRect, etc.)
