
#include <math.h>
/**
 * @brief error function based on the Abramowitz and Stegun approximation
 * 
 * @param x input number
 * @return double 
 */
double myerf01B(double x);
/**
 * @brief inverse error function based on the mclaurin series, with 5 coeficients
 * 
 * @param x input function number
 * @return double 
 */
double myerfinv03B(double x);