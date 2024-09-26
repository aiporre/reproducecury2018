#define Cauchy 1
#define Gauss 2

#if FUN == Cauchy
    #define FUNEVAL(u) 1.0/(1.0+u)
    #define FUNDIFF(u) -1.0/((1.0+u)*(1.0+u))
    #define FUNDIFF2(u) 2.0/((1.0+u)*(1.0+u)*(1.0+u))
#elif FUN == Gauss
    #define FUNEVAL(u) exp(-u)
    #define FUNDIFF(u) -exp(-u)
    #define FUNDIFF2(u) exp(-u)
#elif FUN == Exp
    #define FUNEVAL(u) exp(-sqrt(u))
    #define FUNDIFF(u) u ? -exp(-sqrt(u))/(2*sqrt(u)) : 0
    #define FUNDIFF2(u) u ? exp(-sqrt(u))*(1+1/sqrt(u))/(4*u) : 0
#else
    #warning "FUN is not defined !!"
#endif


