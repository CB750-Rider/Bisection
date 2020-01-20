/* bisection.h is a generic library for performing bisection of
 * a given function. */

include <error.h>

double biscection_solve(double (*y)(double *x, void *args),
		double y0);


