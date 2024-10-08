#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dopri8.h"

#define typ double

typ epsi=0.00000085728263;
typ eta=2.0*M_PI;

void res_346(int n, typ t, typ *y, typ *dy){

	dy[1]=1.5*eta*epsi*y[1]+(1.0+epsi)*eta*y[2]; //valeurs propres %i+epsi et -%i+epsi + O(%i*epsi**2)
	dy[2]=(epsi-1.0)*eta*y[1]+0.5*eta*epsi*y[2];

}

void solout(typ t, const typ *y)
{
	printf("%.3lf   %.13lf   %.13lf\n",t,y[1],y[2]);
}

int main(int argc, char *argv[])
{
	typ tend;
	
	int k=1;
	typ y10=-0.002715, y20=0.00361783;
	
	typ y[3]= {0,y10,y20};  // initial condition of the system
	int n=2;                       // number of dimension of the system
	typ t0=0.0;                // initial time
	typ eps=5.0E-14;              // precision of the integration
	typ hmax=0.1;             // maximal step size
	typ h=0.5*hmax;               // local precision of the integration
	
	int nstep =100000;               // number of steps
	typ stepsize =0.56;//30.0*caractime/((typ) nstep);          // time step for the output
	typ t = t0;

	printf("epsilon=%.13lf,   tps_carac=%.13lf\n\n\n",epsi,1.0/(2*M_PI*epsi));
   
	int dope=0;
	while(k<=nstep && !dope)
    {
     	tend = t0+k*stepsize;
		solout(t, y);
        	dope=gdopri8(&n,trivial,y, &t, &tend, &eps, &hmax, &h);
		k=k+1;
    }
    return 0;    
}
