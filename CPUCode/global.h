#ifndef GLOBAL_H
#define GLOBAL_H

#ifdef TEST_MAIN
#  define EXTERN
#else
#  define EXTERN extern
#endif

EXTERN int T, VOLUME;
EXTERN int LX, LY, LZ;

#endif
