#include "TargetConditionals.h"

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
#pragma clang diagnostic ignored "-Wmissing-prototypes"

#if TARGET_OS_IPHONE // TV, IOS, SIM
#define CR_HIDDEN_SYMBOL_LINK __private_extern__
#else // MAC
#define CR_HIDDEN_SYMBOL_LINK static inline
#endif

CR_HIDDEN_SYMBOL_LINK int __cr_stub_func(void) {
	//#define VARCONCAT(PREFIX, NAME) PREFIX ## NAME
	//	static int VARCONCAT(__VAR, CREACEED_STUB_FUNCNAME) = 0;
	return 0;
}

#pragma clang diagnostic pop
#pragma clang diagnostic pop
