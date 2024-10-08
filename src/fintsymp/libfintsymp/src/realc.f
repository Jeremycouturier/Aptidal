! definition REALC et valeur absolue du treal
#if TREAL == 8 
#define ABSR(x) ABS(x)
#ifdef __GFORTRAN__
#define REALC(x) (x/**/_8)
#else
#define REALC(x) (x##_8)
#endif

#elif TREAL == 10 
#ifdef __GFORTRAN__
#define REALC(x) (x/**/_10)
#else
#define REALC(x) (x##_10)
#endif
#define ABSR(x) ABS(x)

#elif TREAL == 16
#ifdef __GFORTRAN__
#define REALC(x) (x/**/_16)
#else
#define REALC(x) (x##_16)
#endif
#define ABSR(x) ABS(x)

#else 
...
#endif

! format du treal
#if TREAL == 8 
#define FMT_TREAL D23.16
#define SFMT_TREAL "D23.16"
#elif TREAL == 10
#define FMT_TREAL E33.16
#define SFMT_TREAL "E33.16"
#elif TREAL == 16
#define FMT_TREAL E33.16
#define SFMT_TREAL "E33.16"
#else 
...
#endif

