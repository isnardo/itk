#ifndef _RANDOM_H
#define _RANDOM_H

#include <math.h>
#include <stdlib.h>

//Generate random number with uniform distribution between 0 to 1
inline float rand_uniform(void){return ((float)rand()/((float)(RAND_MAX)+(float)(1)));};

//Generate random number with normal distribution for a given variance
float rand_normal(float var){
	float U1,U2,V1,V2,X1;
	float S = 2.0;
	if( var==0 )
		return 0;
	else{
		while( S>=1 ){
			U1 = rand_uniform();
			U2 = rand_uniform();
			V1 = 2.0*U1 - 1.0;
			V2 = 2.0*U2 - 1.0;
			S = V1*V1 + V2*V2;
		}
		X1 = V1*sqrt((-2.0*log(S))/S);
		X1 *= sqrt(var);
		//X2=(1.0/sqrt(2*M_PI*var))*exp(-(X1*X1)/(2*var));
		return X1;
	}
}

#endif
