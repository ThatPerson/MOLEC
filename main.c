#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Infinity 100

typedef float (*Equation)(float*, float);


float test_case(float *args, float x) {
	return args[1];
}

float R(float *args, float r) {
	int n = (int) args[0];
	int l = (int) args[1];
	int Z = (int) args[2];
	//printf("%d %d %d %f\n", n, l, Z, r);
	if (l >= n || l < 0) {
		printf("0 <= l < n\n");
		exit(-1);
	}
	//return exp(-r);
	float p = Z * r;
	if (n == 1) {
		if (l == 0) {
			return 2 * pow(Z, 3./2.) * exp(-p);
		}
	} else if (n == 2) {
		if (l == 1) {
			return pow(Z, 3./2.) * pow(24, -1/2.) * p * exp(-p/2);
		} else if (l == 0) {
			return pow(Z, 3./2.) * pow(8, -1/2.) * (2-p) * exp(-p/2);
		}
	} else if (n == 3) {
		//printf("%f\n",  pow(Z, 3./2.));
		if (l == 2) {
			return pow(Z, 3./2.) * (4/81.) * pow(30, -1/2.) * pow(p, 2) * exp(-p/3);
		} else if (l == 1) {
			return pow(Z, 3./2.) * (4/81.) * pow(6, -1/2.) * (6*p - pow(p, 2)) * exp(-p/3);
		} else if (l == 0) {
			return pow(Z, 3./2.) * (2/81.) * pow(3, -1/2.) * (27 - 18*p + 2*pow(p, 2)) * exp(-p/3);
		}
	} else {
		printf("Don't have this Radial\n");
		exit(-1);
	}
	return 1;
}

float A(float *args, float theta, float phi) {
}
		

float integrate(Equation f, float lower_range, float upper_range, float step, float * args) {
	float l;
	float sum = 0;
	for (l = lower_range; l <= upper_range; l += step) {
		sum += pow(l * f(args, l), 2) * step; 
		//printf("%f\n", f(args, l));
		//printf("%f %f\n", sum, f(args, l));
	}
	return sum;
}

float differentiate(Equation f, float pos, float accuracy, float *args) {
	float upper = pos + (accuracy / 2);
	float lower = pos - (accuracy / 2);
	return (f(args,upper) - f(args, lower))/accuracy;
}

int main(void) {
	float ar[4] = {1, 2, 3, 4};
	//printf("%f\n", test_case(ar, 3));
//	printf("%f\n", integrate(&test_case, 0, 1, 0.0001, ar));
	//printf("%f\n", differentiate(&test_case, 3, 0.0001, ar));

	float n, l;

	for (ar[0] = 0; ar[0] <= 3; ar[0]++) {
		for (ar[1] = 0; ar[1] < ar[0]; ar[1]++) {
			
			ar[2] = 1;
			printf("%f %f: %f\n", ar[0], ar[1], integrate(R, 0, Infinity, 0.1, ar));
		}
	}
	return 1;
}
