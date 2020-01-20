/*Copyright (c) 2020 Lee R. Burchett

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/

#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <errno.h>
#include <float.h>
#include "bisection.h"

#define MAX_ITERATIONS 1000


static bool get_next_x(double *x);
static bool get_start(
	double (*y)(double x, void *args),
	double y0,
	double x[N_A],
	void *args);
static double run_bisection_iteration(
	double (*y)(double x, void *args),
	double y0,
	double x[N_A],
	void *args);
static bool is_almost(double x, double y);
static bool check_if_converged(
		double y0,
    double y,
		double x[N_A]);

/* Bisection_solve is a public entry-point function
 * for using the bisection method. Bisection_solve takes two
 * required arguments (y ,nd y0) and two optional arguments
 * (args, x_init). It returns the value of x that gives y0.
 *
 * Inputs:
 * y      a function that takes in a value of x and returns a y value.
 * y0     the value for which we wish to find a corresponding x0.
 * args   a pointer to the function arguments.
 * x_init an optional initial x range (ABOVE, BELOW, CENTER). 
 *        pass NULL if you wish for x_init to be determined automatically. 
 *
 * Outputs:
 * x0     the x value that gives y0. */
double bisection_solve(
		double (*y)(double x, void *args),
		double y0,
    void *args,
    double *x_init)
{
  double x[N_A]={0.0},yc=0.0;
  int watchdog=0;

  if(x_init==NULL){
    if(get_start(y,y0,x,args)){
      errno=EINVAL;
      return 0.0;
    }
  } /* if x_init is NULL*/
  else {
    x[ABOVE] = x_init[ABOVE];
    x[BELOW] = x_init[BELOW];
    x[CENTER] = x_init[CENTER];
  } /* else copy the start values */

  do{
  	yc = run_bisection_iteration(y,y0,x,args);
  	if(watchdog++>MAX_ITERATIONS){
  		errno=EINVAL;
  		return 0.0;
  	}
  } while(!check_if_converged(yc,y0,x));

  return x[CENTER];
} /* bisection_solve */

/* get_next_x
 * TODO: This function needs some thought to
 * get the implimentation right. What is a good
 * way to explore the available space of doubles?
 * For now, create a Fibbinacci sequence. */
bool get_next_x(double *x){
  static double next_x=1.0;
  static double last_x=0.0;
  if(next_x>0.0){
    if( (FLT_MAX-next_x) <= last_x){ /* We will overflow next time. */
      *x = next_x;
      next_x = -1.0; /* Start running negative values. */
      last_x = 0.0;
    } /* if we will overflow next time */
    else{
      last_x = *x;
      *x = next_x;
      next_x += last_x;
    } /* else we won't overflow  */
  }/* if next_x > 0 */ 
  else{
    if( (FLT_MAX+next_x) <= -last_x){ /* We will overflow next time. */
      *x = next_x;
      /* End of the line */
      return true;
    } /* if we will overflow next time */
    else{
      last_x = *x;
      *x = next_x;
      next_x += last_x;
    } /* else we won't overflow  */
  }/* else next_x < 0.0*/
  return false;
}/* get_next_x*/


/* Get_start attempts to find a starting set of values for x. */
bool get_start(
	double (*y)(double x, void *args),
	double y0,
	double x[N_A],
	void *args)
{
  x[BELOW] = 0.0;
  while(y(x[BELOW],args) > y0){
    if(get_next_x(&x[BELOW])){
      /* get_next_x only returns true if we hit our watchdog condition. 
       * This indicates an error as no value of x was found that
       * gives us a y value less than y0. */
      return true; /* Alert to an error. */
    } /* if we get a watchdog error. */
  } /* while y(x[BELOW]) > y0 */
  x[ABOVE] = 0.0;
  while(y(x[ABOVE],args) < y0){
    if(get_next_x(&x[ABOVE])){
      return true;
    } /* if we get a watchdog error */
  } /* while y(x[ABOVE]) < y0 */
  x[CENTER] = 0.5*(x[BELOW] + x[ABOVE]);
  return false;
} /* get_start */

/* Run_bisection_iteration performs one iteration of the bisection 
 * method. It updates the x estimates and returns the most recent
 * value of y(x[CENTER]). */
double run_bisection_iteration(
	double (*y)(double x, void *args),
	double y0,
	double x[N_A],
	void *args)
{
  /* A new yc is calculated as the function value at x[CENTER]. */
	double yc = y(x[CENTER],args);
  /* Figure out if yc is > or less than our target y0.
   * Change either X[ABOVE] or X[BELOW] to be equal to the current
   * x[CENTER] and so that y(x[ABOVE]) >= y0 >= y(x[BELOW]). */
	if(yc > y0){
		x[ABOVE] = x[CENTER];
	}
	else{
		x[BELOW] = x[CENTER];
	}
	x[CENTER] = 0.5*(x[ABOVE] + x[BELOW]); /* x[CENTER] for the 
  next iteration is half-way betweer x[ABOVE] and x[BELOW]. */
	return yc;
} /* run_bisection_iteration */

bool is_almost(double x, double y)
{
  /* Check to see if two numbers are nearly identical. */
  double diff = fabs(x-y); 

  /* Handle some casse that could give us trouble. */
  if(diff==0.0) return true;
  if(x==0.0) return fabs(y) < DBL_MIN*2.0;
  if(y==0.0) return fabs(x) < DBL_MIN*2.0;

  /* Convert everything to log2 */
  diff = log2(diff);
  x = log2(fabs(x));
  y = log2(fabs(y));

  if (x > y){
    return ((y-diff) > (DBL_MANT_DIG - 2)); /* Is our difference less than machine precision for both x and y? */
  } /* if x > y */
  else {
    return ((x-diff) > (DBL_MANT_DIG - 2)); 
  } /* else y > x */
} /* is_almost */

bool check_if_converged(
		double y0,
		double yc,
		double x[N_A])
{
  /* Check to see if the difference in y or x is now so
  * small that we can consider things as "converged". */
  if(is_almost(y0,yc) || is_almost(x[ABOVE],x[BELOW]))
	  return true;
  else
	  return false;
} /* check_if_converged */
