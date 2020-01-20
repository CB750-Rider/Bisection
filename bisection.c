#define MAX_ITERATIONS 1000

typedef struct {
	Abv,
	Blw,
	C,
	N_A} ESTIMATES;

int has_not_converged(
		double (*y)(double *x,void *args),
		double y0,
		double x[N_A]);

double biscection_solve(
		double (*y)(double *x, void *args),
		double y0)
{
double x[N_A],yc;
int watchdog=0;

get_start(y,y0,x);

do{
	yc = bisect(y,y0,x);
	if(watchdog++>MAX_ITERATIONS){
		set errno=EINVAL;
		return 0.0;
	}

} while(has_not_converged(y,y0,x));

return x[C];
}

double bisect(
	double (*y)(double x, void *args),
	double y0,
	double x[N_A],
	void *args)
{
	double yc = y(x[C],args);
	if(yc > y0){
		x[Abv] = x[C];
	}
	else{
		x[Blw] = x[C];
	}
	x[C] = 0.5*(x[Abv] + x[Blw]);
	return yc;
}
int has_not_converged(
		double y0,
		double yc,
		double x[N_A])
{
if(check_for_convergence(y0,yc) ||
  check_for_convergence(x[Abv],x[Blw]))
	return 0;
else
	return 1;
}
