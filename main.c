#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>

#define Infinity 100
#define PI 3.14159265358

#define hbar 1.05457168 * pow(10, -34)
#define Me 0.91093826 * pow(10, -30)
#define ec 1.60217653 * pow(10, -19) 
#define e0 8.8541878 * pow(10, -12)


typedef float (*Equation)(float*, float);
typedef float (*Wavefunction)(int, int, int, float, float, float, float, float, float);

float det2(float H[2][2], float S[2][2], float E) {
	//[ H00 - ES00	H10 - ES10 ]
	//[ H01 - ES01	H11 - ES11 ]
	
	return ((H[1][1] - E*S[1][1]) * (H[0][0] - E*S[0][0])) - ((H[1][0] - E*S[1][0]) * (H[0][1] - E*S[0][1]));
}

float det3(float H[3][3], float S[3][3], float E) {
	float fe = (H[0][0] - E*S[0][0]) * (((H[1][1] - E*S[1][1]) * (H[2][2] - E*S[2][2])) - ((H[2][1] - E*S[2][1]) * (H[1][2] - E*S[1][2])));
	float fi = (H[1][0] - E*S[1][0]) * (((H[0][1] - E*S[0][1]) * (H[2][2] - E*S[2][2])) - ((H[2][1] - E*S[2][1]) * (H[0][2] - E*S[0][2])));
	float fo = (H[2][0] - E*S[2][0]) * (((H[0][1] - E*S[0][1]) * (H[1][2] - E*S[1][2])) - ((H[0][2] - E*S[0][2]) * (H[1][1] - E*S[1][1])));
	float fum = fe + fi + fo;
	return fum;
}	


float wavefunction(int n, int l, int q, float pos[3], float offset[3], int Z) {
	// args is [n, l, m, Z]
	float x = pos[0] - offset[0];
	float y = pos[1] - offset[1];
	float z = pos[2] - offset[2];
	float r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
	if (n == 1) {
		return pow(Z, 3/2.) * pow(1/PI, 1/2.) * exp(-r);
//	} else if (strcmp(name, "2s") == 0) {
	} else if (n == 2) {
		if (l == 0) {
			return pow(Z, 3/2.) * (1/2)*pow(1/(2*PI), 1/2.) * (1 - r/2) * exp(-r/2);
		} else if (l == 1) {
			if (q == 0) {
				return pow(Z, 3/2.) * (1 / (4*sqrt(2*PI))) * z * exp(-r/2);
			} else if (q == 1) {
				return pow(Z, 3/2.) * (1 / (4*sqrt(2*PI))) * x * exp(-r/2);
			} else if (q == 2) {
				return pow(Z, 3/2.) * (1 / (4*sqrt(2*PI))) * y * exp(-r/2);
			}
		}
	} else if (n == 3) {
		if (l == 0) {
			return pow(Z, 3/2.) * (1 / (81*sqrt(3*PI))) * (27 - 18*r + 2*pow(r, 2)) * exp(-r/3);
		} else if (l == 1) {
			if (q == 0) {
				return pow(Z, 3/2.) * (sqrt(2) / (81*sqrt(PI))) * (6 - r) * z * exp(-r/3);
			} else if (q == 1) {
				return pow(Z, 3/2.) * (sqrt(2) / (81*sqrt(PI))) * (6 - r) * y * exp(-r/3);
			} else if (q == 2) {
				return pow(Z, 3/2.) * (sqrt(2) / (81*sqrt(PI))) * (6 - r) * x * exp(-r/3);
			}
		} else if (l == 2) {
			if (q == 0) {
				return pow(Z, 3/2.) * (1 / (81*sqrt(6 * PI))) * (3 * pow(z, 2) - pow(r, 2)) * exp(-r/3);
			} else if (q == 1) {
				return pow(Z, 3/2.) * (sqrt(2) / (81 * sqrt(PI))) * z * x * exp(-r/3);
			} else if (q == 2) {
				return pow(Z, 3/2.) * (sqrt(2) / (81 * sqrt(PI))) * z * y * exp(-r/3);
			} else if (q == 3) {
				return pow(Z, 3/2.) * (sqrt(2) / (81 * sqrt(PI))) * x * y * exp(-r/3);
			} else if (q == 4) {
				return pow(Z, 3/2.) * (1 / (81 * sqrt(PI))) * (pow(x, 2) - pow(y, 2)) * exp(-r/3);
			}
		}
	} else {
		printf("I don't recognise that orbital..\n");
		exit(1);
	}
}

float Hamiltonian(int nA, int lA, int qA, float offsetA[3], float posA[3], int Z) {
	// Numerical Form of Hamiltonian
	float alpha = 0.0001;
	float beta = 0.00001; // Beta needs to be smaller than alpha
	// r component
	
	float r = sqrt(pow(posA[0], 2) + pow(posA[1], 2) + pow(posA[2], 2));
	
	float dxsquu[4] = {posA[0] + alpha + beta, posA[1], posA[2]};
	float dxsqud[4] = {posA[0] + alpha - beta, posA[1], posA[2]};
	float dxsqdu[4] = {posA[0] - alpha + beta, posA[1], posA[2]};
	float dxsqdd[4] = {posA[0] - alpha - beta, posA[1], posA[2]};
	
	float dysquu[4] = {posA[0], posA[1] + alpha + beta, posA[2]};
	float dysqud[4] = {posA[0], posA[1] + alpha - beta, posA[2]};
	float dysqdu[4] = {posA[0], posA[1] - alpha + beta, posA[2]};
	float dysqdd[4] = {posA[0], posA[1] - alpha - beta, posA[2]};
	
	float dzsquu[4] = {posA[0], posA[1], posA[2] + alpha + beta};
	float dzsqud[4] = {posA[0], posA[1], posA[2] + alpha - beta};
	float dzsqdu[4] = {posA[0], posA[1], posA[2] - alpha + beta};
	float dzsqdd[4] = {posA[0], posA[1], posA[2] - alpha - beta};
	
	float dxsq = (((wavefunction(nA, lA, qA, dxsquu, offsetA, Z) - wavefunction(nA, lA, qA, dxsqud, offsetA, Z)) / (2*beta)) - ((wavefunction(nA, lA, qA, dxsqdu, offsetA, Z) - wavefunction(nA, lA, qA, dxsqdd, offsetA, Z)) / (2*beta))) / (2*alpha);
	
	float dysq = (((wavefunction(nA, lA, qA, dysquu, offsetA, Z) - wavefunction(nA, lA, qA, dysqud, offsetA, Z)) / (2*beta)) - ((wavefunction(nA, lA, qA, dysqdu, offsetA, Z) - wavefunction(nA, lA, qA, dysqdd, offsetA, Z)) / (2*beta))) / (2*alpha);
		
	float dzsq = (((wavefunction(nA, lA, qA, dzsquu, offsetA, Z) - wavefunction(nA, lA, qA, dzsqud, offsetA, Z)) / (2*beta)) - ((wavefunction(nA, lA, qA, dzsqdu, offsetA, Z) - wavefunction(nA, lA, qA, dzsqdd, offsetA, Z)) / (2*beta))) / (2*alpha);
	
	float kineticenergy = - (dxsq + dysq + dzsq);
	
	float potentialenergy = - Z / r;
	
	float totalenergy = kineticenergy + potentialenergy;
	

	
	//printf("%f %f %f %f\n", rse, thetase, nA, lA, qAse, chase);
	
	return totalenergy;
}

float calculate_energy(int nA, int lA, int qA, int nB, int lB, int qB, int Za, int Zb, float bond_length, int siz) {
	float numerator = 0;
	float denominator = 0;

	float hammy = 0;
	float AWav = 0;
	float BWav = 0;
	
	float x = 0;
	float y = 0; 
	float z = 0;
	float sA[3] = {bond_length, 0, 0};
	float sB[3] = {0, 0, 0};
	//float pos[3];
	float step = 0.1;
	
	/*int i, thread_id;
    int glob_nloops, priv_nloops;
    glob_nloops = 0;

    // parallelize this chunk of code
    #pragma omp parallel private(priv_nloops, thread_id) reduction(+:glob_nloops)
    {
        priv_nloops = 0;
        thread_id = omp_get_thread_num();

        // parallelize this for loop
        #pragma omp for
        for (i=0; i<100000; ++i)
        {
            ++priv_nloops;
        }
        glob_nloops += priv_nloops;
    }
    printf("Total # loop iterations is %d\n",
           glob_nloops);*/
	int pos_x, pos_y, pos_z;
	float pos[3];
	int priv_nloops, priv_num, priv_den;
	int thread_id;
	#pragma omp parallel private(priv_num, priv_den, pos) reduction(+:numerator,denominator)
	{
		priv_nloops = 0;
		priv_num = 0;
		priv_den = 0;
		#pragma omp for
		for (pos_x = 0; pos_x <= 1000; pos_x ++) {
			for (pos_y = 0; pos_y <= 1000; pos_y ++) {
				for (pos_z = 0; pos_z <= 1000; pos_z ++) {
					//printf("%f %f %f\n", pos[0], pos[1], pos[2]);
					pos[0] = (pos_x - 500)/10.;
					pos[1] = (pos_y - 500)/10.;
					pos[2] = (pos_z - 500)/10.;
					thread_id = omp_get_thread_num();
					hammy = Hamiltonian(nA, lA, qA, sA, pos, Za);
					AWav = wavefunction(nA, lA, qA, pos, sA, Za);
					BWav = wavefunction(nB, lB, qB, pos, sB, Zb);
					//printf("%f %f %f done by thread %d\n", pos[0], pos[1], pos[2], thread_id);
					if (hammy != NAN && AWav != NAN && BWav != NAN) {

						priv_num += hammy * BWav * pow(0.1, 3);
						priv_den += AWav * BWav * pow(0.1, 3);
					} else {
						printf("It's non existent?\n");
					}
				}
			}
		}
		numerator += priv_num;
		denominator += priv_den;
	}
	float energy = numerator/denominator;
	//printf("%f %f\n", numerator, denominator);
	return energy;
	
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



	/*for (ar[0] = 0; ar[0] <= 3; ar[0]++) {
		for (ar[1] = 0; ar[1] < ar[0]; ar[1]++) {
			
			ar[2] = 1;
			printf("%f %f: %f\n", ar[0], ar[1], integrate(R, 0, Infinity, 0.1, ar));
		}
	}*/
	
	float rs[4] = {2, 0, 0, 1};
	int n = 50;
// float calculate_energy(char * nA, lA, qA, char * nB, lB, qB, int Za, int Zb, float bond_length) {
	printf("%d, %f\n", n, calculate_energy(1, 0, 0, 1, 0, 0, 1, 1, 0, n));
	printf("%d, %f\n", n, calculate_energy(2, 0, 0, 2, 0, 0, 1, 1, 0, n));
	return 1;
}
