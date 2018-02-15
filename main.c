#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>

#define Infinity 100
#define PI 3.14159265358
#define d_size 50
#define d_step 1
#define hbar 1.05457168 * pow(10, -34)
#define Me 0.91093826 * pow(10, -30)
#define ec 1.60217653 * pow(10, -19) 
#define e0 8.8541878 * pow(10, -12)

typedef struct {
	double a;
	double b;
} pair;
typedef double (*Equation)(double*, double);
typedef double (*Wavefunction)(int, int, int, double, double, double, double, double, double);

float det(float m[50][50], int dimensions);
float det(float m[50][50], int dimensions) {
	int i, x, y, a, b;
	if (dimensions == 1) {
		return m[0][0];
	}
	
	float o[50][50];
	float prod = 0;
	
	int max_dim = dimensions - 1;
	for (i = 0; i < dimensions; i++) {
		a = 0; b = 0;
		for (x = 1; x < dimensions; x++) {	
			for (y = 0; y < dimensions; y++) {
				if (y != i) {
					o[a][b] = m[x][y];
					
					a++;
					if (a > max_dim-1) {
						a = 0;
						b++;
					}
					
				}
				
			}
		}
		prod = prod + pow(-1, i) * m[0][i] * det(o, dimensions - 1);
	}
	return prod;
}

double det2(double H[2][2], double S[2][2], double E) {
	//[ H00 - ES00	H10 - ES10 ]
	//[ H01 - ES01	H11 - ES11 ]
	
	return ((H[1][1] - E*S[1][1]) * (H[0][0] - E*S[0][0])) - ((H[1][0] - E*S[1][0]) * (H[0][1] - E*S[0][1]));
}

double det3(double H[3][3], double S[3][3], double E) {
	double fe = (H[0][0] - E*S[0][0]) * (((H[1][1] - E*S[1][1]) * (H[2][2] - E*S[2][2])) - ((H[2][1] - E*S[2][1]) * (H[1][2] - E*S[1][2])));
	double fi = (H[1][0] - E*S[1][0]) * (((H[0][1] - E*S[0][1]) * (H[2][2] - E*S[2][2])) - ((H[2][1] - E*S[2][1]) * (H[0][2] - E*S[0][2])));
	double fo = (H[2][0] - E*S[2][0]) * (((H[0][1] - E*S[0][1]) * (H[1][2] - E*S[1][2])) - ((H[0][2] - E*S[0][2]) * (H[1][1] - E*S[1][1])));
	double fum = fe + fi + fo;
	return fum;
}	



double wavefunction(int n, int l, int q, double pos[3], double offset[3], int Z) {
	// args is [n, l, m, Z]
	double x = pos[0] - offset[0];
	double y = pos[1] - offset[1];
	double z = pos[2] - offset[2];
	//printf("%lf %lf %lf\n", x, y, z);
	double r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
	if (n == 1) {
		return pow(Z, 3/2.) * pow(1/PI, 1/2.) * exp(-r);
//	} else if (strcmp(name, "2s") == 0) {
	} else if (n == 2) {
		if (l == 0) {
			//printf("%lf %lf %lf %lf == %lf\n", pow(Z, 3/2.), (1/2.)*pow(1/(2*PI), 1/2.) ,  (1 - r/2), exp(-r/2), pow((double) Z, 3/2.) * (1/2.)*pow(1/(2*PI), 1/2.) * (1 - r/2) * exp(-r/2));
			return pow((double) Z, 3/2.) * (1/2.)*pow(1/(2*PI), 1/2.) * (1 - r/2) * exp(-r/2);
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

double Hamiltonian(int nA, int lA, int qA, double offsetA[3], double posA[3], int Z) {
/* This is the problem!!! */


	// Numerical Form of Hamiltonian
	double alpha = 0.01;
	double beta = 0.001; // Beta needs to be smaller than alpha
	// r component
	//posA[0] = posA[0]/10;
	//posA[1] = posA[1]/10;
	//posA[2] = posA[2]/10;
	double r = sqrt(pow(posA[0], 2) + pow(posA[1], 2) + pow(posA[2], 2));
	if (r == 0) {
		// We are at the center, and so the potential will make it fail.
		return NAN;
	}
	
	double dxsquu[4] = {posA[0] + alpha + beta, posA[1], posA[2]};
	double dxsqud[4] = {posA[0] + alpha - beta, posA[1], posA[2]};
	double dxsqdu[4] = {posA[0] - alpha + beta, posA[1], posA[2]};
	double dxsqdd[4] = {posA[0] - alpha - beta, posA[1], posA[2]};
	
	double dysquu[4] = {posA[0], posA[1] + alpha + beta, posA[2]};
	double dysqud[4] = {posA[0], posA[1] + alpha - beta, posA[2]};
	double dysqdu[4] = {posA[0], posA[1] - alpha + beta, posA[2]};
	double dysqdd[4] = {posA[0], posA[1] - alpha - beta, posA[2]};
	
	double dzsquu[4] = {posA[0], posA[1], posA[2] + alpha + beta};
	double dzsqud[4] = {posA[0], posA[1], posA[2] + alpha - beta};
	double dzsqdu[4] = {posA[0], posA[1], posA[2] - alpha + beta};
	double dzsqdd[4] = {posA[0], posA[1], posA[2] - alpha - beta};
	
	//		dxsquu[0], dxsquu[1], dxsquu[2], dxsquu[3], 
	//		dxsqud[0], dxsqud[1], dxsqud[2], dxsqud[3], 
	//		dxsqdu[0], dxsqdu[1], dxsqdu[2], dxsqdu[3], 
	//		dxsqdd[0], dxsqdd[1], dxsqdd[2], dxsqdd[3]);
	
	// -(2.64589 Power[10, -10]) (D[D[Hx[x, y, z], x], x] - 
   //  D[D[Hx[x, y, z], y], y] - D[D[Hx[x, y, z], z], z]) - 
   //  10/Sqrt[x^2 + y^2 + z^2] Hx[x, y, z];
	
	
	
	/* Work out how to map this onto the actual values.!!!!!!!!!!!!!! */
	
	
	double dxsq = (((wavefunction(nA, lA, qA, dxsquu, offsetA, Z) - wavefunction(nA, lA, qA, dxsqud, offsetA, Z)) / (2*beta)) - ((wavefunction(nA, lA, qA, dxsqdu, offsetA, Z) - wavefunction(nA, lA, qA, dxsqdd, offsetA, Z)) / (2*beta))) / (2*alpha);
	double dysq = (((wavefunction(nA, lA, qA, dysquu, offsetA, Z) - wavefunction(nA, lA, qA, dysqud, offsetA, Z)) / (2*beta)) - ((wavefunction(nA, lA, qA, dysqdu, offsetA, Z) - wavefunction(nA, lA, qA, dysqdd, offsetA, Z)) / (2*beta))) / (2*alpha);
	double dzsq = (((wavefunction(nA, lA, qA, dzsquu, offsetA, Z) - wavefunction(nA, lA, qA, dzsqud, offsetA, Z)) / (2*beta)) - ((wavefunction(nA, lA, qA, dzsqdu, offsetA, Z) - wavefunction(nA, lA, qA, dzsqdd, offsetA, Z)) / (2*beta))) / (2*alpha);
	
	//printf("%lf %lf %lf :SQUARES: %lf %lf %lf\n", posA[0], posA[1], posA[2], dxsq, dysq, dzsq);
	
	double a = 6.10426 * pow(10, 1);
	double b = 2.30708 * pow(10, 12);
	
	double kineticenergy = - a * (dxsq + dysq + dzsq);
	double potentialenergy = - b * wavefunction(nA, lA, qA, posA, offsetA, Z) * Z / r;
	double totalenergy = kineticenergy + potentialenergy;
	
	
	
	//printf("%lf %lf %lf %lf\n", rse, thetase, nA, lA, qAse, chase);
	
	return (double)totalenergy * (2.1673 * pow(10, -13));
}

pair calculate_energy(int nA, int lA, int qA, int nB, int lB, int qB, int Za, int Zb, double bond_length, int siz) {

	/* My current assumption would be that you'd need two functions - H, and S. H does the Hermitian part, S does the overlap part. Then you can directly use those in the MO calculations, and in normal use can just to like Hii/Sii */

	double numerator = 0;
	double denominator = 0;

	double hammy = 0;
	double AWav = 0;
	double BWav = 0;
	
	double x = 0;
	double y = 0; 
	double z = 0;
	double sA[3] = {bond_length, 0, 0};
	double sB[3] = {0, 0, 0};
	double step = 0.1;
	double pos_x, pos_y, pos_z;
	double pos[3];
	int priv_nloops;
	double priv_num, priv_den;
	int thread_id;
	//#pragma omp parallel private(priv_num, priv_den, pos) reduction(+:numerator,denominator)
	//{
		priv_nloops = 0;
		priv_num = 0;
		priv_den = 0;
		double r, theta, phi;
		//#pragma omp for
		for (r = 0.01; r <= 40; r+=0.01) {
			for (theta = 0.01; theta <= PI; theta+=0.1) {
				for (phi = 0.01; phi <= 2 * PI; phi+=0.1) {
					pos[0] = r * cos(phi) * sin(theta);
					pos[1] = r * sin(phi) * sin(theta);
					pos[2] = r * cos(theta);
					hammy = Hamiltonian(nA, lA, qA, sA, pos, Za);
					AWav = wavefunction(nA, lA, qA, pos, sA, Za);
					BWav = wavefunction(nB, lB, qB, pos, sB, Zb);
					if (!isnan(hammy) && !isnan(AWav) && !isnan(BWav)) {
						priv_num += hammy * BWav * pow(r, 2) * sin(theta);
						priv_den += AWav * BWav * pow(r, 2) * sin(theta); 
					} else {
						printf("It's non existent?\n");
					}
				}
			}
		}
		numerator += priv_num;
		denominator += priv_den;
	//}
	printf("%lf %lf\n", numerator, denominator);
	pair xl;
	xl.a = numerator;
	xl.b = denominator;
	double energy = numerator/denominator;
	//printf("%lf %lf\n", numerator, denominator);
	return xl;
	
}



double integrate(Equation f, double lower_range, double upper_range, double step, double * args) {
	double l;
	double sum = 0;
	for (l = lower_range; l <= upper_range; l += step) {
		sum += pow(l * f(args, l), 2) * step; 
		//printf("%lf\n", f(args, l));
		//printf("%lf %lf\n", sum, f(args, l));
	}
	return sum;
}

double differentiate(Equation f, double pos, double accuracy, double *args) {
	double upper = pos + (accuracy / 2);
	double lower = pos - (accuracy / 2);
	return (f(args,upper) - f(args, lower))/accuracy;
}

int main(void) {
	double ar[4] = {1, 2, 3, 4};
	double rs[4] = {2, 0, 0, 1};
	int n = 50;
	
	
	/*double hammy, r, theta, phi;		
	double pos[3];
	double s[3] = {0, 0, 0};
	for (theta = 0; theta < 3.14; theta += 0.1) {
		r = 2;
		//theta = 1;
		phi = 1;
		pos[0] = r * cos(phi) * sin(theta);
		pos[1] = r * sin(phi) * sin(theta);
		pos[2] = r * cos(theta);
		hammy = Hamiltonian(1, 0, 0, s, pos, 1);
		printf("%lf, %lf\n", theta, hammy);
	}*/
	
	long double ns = 0.000000001;
	printf("%Lf\n", ns);
	
	int orbitals[3][3];
	orbitals[0][0] = 1;
	orbitals[0][1] = 0;
	orbitals[0][2] = 0;
	orbitals[1][0] = 1;
	orbitals[1][1] = 0;
	orbitals[1][2] = 0;
	orbitals[2][0] = 1;
	orbitals[2][1] = 0;
	orbitals[2][2] = 0;
	
	
	/*printf("1s 1s %d, %lf\n", n, calculate_energy(1, 0, 0, 1, 0, 0, 1, 1, 0, n));
	printf("2s 2s %d, %lf\n", n, calculate_energy(2, 0, 0, 2, 0, 0, 1, 1, 0, n));
	printf("3s 3s %d, %lf\n", n, calculate_energy(3, 0, 0, 3, 0, 0, 1, 1, 0, n));
	printf("2p 2p %d, %lf\n", n, calculate_energy(2, 1, 0, 2, 1, 1, 1, 1, 0, n));
	printf("3p 3p %d, %lf\n", n, calculate_energy(3, 1, 0, 3, 1, 0, 1, 1, 0, n));
	printf("3d 3d %d, %lf\n", n, calculate_energy(3, 2, 0, 3, 2, 0, 1, 1, 0, n));
	
	printf("1s 1s BL1.6 %d, %lf\n", n, calculate_energy(1, 0, 0, 1, 0, 0, 1, 1, 1.6 * pow(10, -8), n));
	printf("1s 1s BL1.4 %d, %lf\n", n, calculate_energy(1, 0, 0, 1, 0, 0, 1, 1, 1.4 * pow(10, -8), n));
	printf("1s 1s BL1.2 %d, %lf\n", n, calculate_energy(1, 0, 0, 1, 0, 0, 1, 1, 1.2 * pow(10, -8), n));
	printf("1s 1s BL1.0 %d, %lf\n", n, calculate_energy(1, 0, 0, 1, 0, 0, 1, 1, 1.0 * pow(10, -8), n));
	printf("1s 1s BL0.8 %d, %lf\n", n, calculate_energy(1, 0, 0, 1, 0, 0, 1, 1, 0.8 * pow(10, -8), n));
	printf("1s 1s BL0.6 %d, %lf\n", n, calculate_energy(1, 0, 0, 1, 0, 0, 1, 1, 0.6 * pow(10, -8), n));
	printf("1s 1s BL0.4 %d, %lf\n", n, calculate_energy(1, 0, 0, 1, 0, 0, 1, 1, 0.4 * pow(10, -8), n));*/
	
	float Ss[5][5];
	float Hs[5][5];
	pair xs;
	int x, y;
	for (x = 0; x < 3; x++) {
		for (y = 0; y < 3; y++) {

			xs = calculate_energy(orbitals[x][0], orbitals[x][1], orbitals[x][2], orbitals[y][0], orbitals[y][1], orbitals[y][2], 1, 1, abs(x - y), n);
			Ss[x][y] = xs.b/10000;
			Hs[x][y] = xs.a/10000;
			
		}
	}
	
	printf("H\n");
	for (x = 0; x < 3; x++) {
		for (y = 0; y < 3; y++) {
			printf("%0.4f\t", Hs[x][y]);
		}
		printf("\n");
	}
	printf("S\n");
	for (x = 0; x < 3; x++) {
		for (y = 0; y < 3; y++) {
			printf("%0.4f\t", Ss[x][y]);
		}
		printf("\n");
	}
	
	
	return 1;
}
