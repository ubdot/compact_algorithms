#include "erf.h"

double myerf01B(double x) {
    double t, rer;
    double ttt;
    double t1 = fabs(x); // Cambiado a fabs para double

     t = 0.078108 * t1 + 0.000972;
    //t = fma(0.078108, t1, 0.000972);

     t = t * t1 + 0.230389;
    //t = fma(t, t1, 0.230389);

     t = t * t1 + 0.278393;
    //t = fma(t, t1, 0.278393);

     t = t * t1 + 1.0;
    //t = fma(t, t1, 1.0);

    ttt = t * t;
    rer = 1.0 - 1.0 / (ttt * ttt);
    
    rer = copysign(rer, x); 

    return rer;
}

double myerfinv03B(double x){
    double temp, xx;
    xx = x * x;

     temp = 0.041147395 * xx + 0.048336063;
    //temp = fma(0.041147395, xx, 0.048336063);

     temp = temp * xx + 0.058372501;
    //temp = fma(temp, xx, 0.058372501);

     temp = temp * xx + 0.073299079;
    //temp = fma(temp, xx, 0.073299079);

     temp = temp * xx + 0.097663620;
    //temp = fma(temp, xx, 0.097663620);

     temp = temp * xx + 0.143931731;
    //temp = fma(temp, xx, 0.143931731);

    temp = temp * xx + 0.261799388;
    //temp = fma(temp, xx, 0.261799388);

    temp = temp * xx + 1.0;
    // temp = fma(temp, xx, 1.0);

    return temp * x * 0.8862;
}