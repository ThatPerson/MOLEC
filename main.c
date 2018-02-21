#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include <glpk.h>


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

typedef struct {
	int orbitals[50][3];
	int num_orbitals;
	float pos[3];
	float zeff;
} atom;

typedef double (*Equation)(double*, double);
typedef double (*Wavefunction)(int, int, int, double, double, double, double, double, double);

float det(float m[50][50], int dimensions);
float det(float m[50][50], int dimensions) {

/* TODO: Calculate minimum energy orbital (H/S for the orbital and itself, for all orbitals present) and maximum. Add half the difference to the top and bottom, and then go over in steps of ~0.01 and look for changes of sign. Mark down and output.
	Also find out what on earth the units are - run the program for a range of distances and then calculate actual for like 100pm in Mathematica and find where it intersects */

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




double wavefunction(int n, int l, int q, double pos[3], double offset[3], float Z) {
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
	return -1;
}

double hamiltonian(int nA, int lA, int qA, double offsetA[3], double posA[3], float Z, float ZB, double posB[3]) {
/* This is the problem!!! */


	// Numerical Form of hamiltonian
	double alpha = 0.01;
	double beta = 0.001; // Beta needs to be smaller than alpha
	// r component
	//posA[0] = posA[0]/10;
	//posA[1] = posA[1]/10;
	//posA[2] = posA[2]/10;
	double r = sqrt(pow(posA[0], 2) + pow(posA[1], 2) + pow(posA[2], 2));
	double rB = sqrt(pow(posB[0] - posA[0], 2) + pow(posB[1] - posA[1], 2) + pow(posB[2] - posA[2], 2));
	if (r == 0) {
		// We are at the center, and so the potential will make it fail.
		return NAN;
	}
	if (rB == 0) {
		rB = 0.001;
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
	
	double a = 1/2.;
	double b = 1;
	
	//double a = 1;
	//double b = 1;
	double kineticenergy = - a * (dxsq + dysq + dzsq);
	double potentialenergy = - b * wavefunction(nA, lA, qA, posA, offsetA, Z) * Z / r;
	double otheratom = - b * wavefunction(nA, lA, qA, posA, offsetA, Z) * ZB / rB;
	//printf("%f %f %f\n", kineticenergy, potentialenergy, otheratom);
	double totalenergy = kineticenergy + (potentialenergy + otheratom);
	
	
	
	//printf("%lf %lf %lf %lf\n", rse, thetase, nA, lA, qAse, chase);
	
	return (double)totalenergy/* * (2.1673 * pow(10, -13))*/;
}

pair calculate_energy(int nA, int lA, int qA, int nB, int lB, int qB, float Za, float Zb, double bond_length, int siz) {

	double numerator = 0;
	double denominator = 0;

	double hammy = 0;
	double AWav = 0;
	double BWav = 0;
	
	double sA[3] = {bond_length, 0, 0};
	double sB[3] = {0, 0, 0};

	double pos[3];

	double priv_num, priv_den;

	//#pragma omp parallel private(priv_num, priv_den, pos) reduction(+:numerator,denominator)
	//{
		priv_num = 0;
		priv_den = 0;
		double r, theta, phi;
		//#pragma omp for
		for (r = 0.01; r <= 40; r+=0.1) {
			for (theta = 0.01; theta <= PI; theta+=0.1) {
				for (phi = 0.01; phi <= 2 * PI; phi+=0.1) {
					pos[0] = r * cos(phi) * sin(theta);
					pos[1] = r * sin(phi) * sin(theta);
					pos[2] = r * cos(theta);
					
			//		double hamiltonian(int nA, int lA, int qA, double offsetA[3], double posA[3], int Z, int ZB, double posB[3]) {
			//		double wavefunction(int n, int l, int q, double pos[3], double offset[3], int Z) {
					
					hammy = hamiltonian(nA, lA, qA, sA, pos, Za, Zb, sB);
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

int switched(float a, float b) {
	if (a >= 0 && b < 0)
		return 1;
	if (b >= 0 && a < 0)
		return 1;
	return 0;
}

void calculate_coefficients(float matr[50], int len, float * arr) {
	/*arr[0] = 3;
	printf("======START======\n");
	f
	// Minimise sum of matr[n] subject to sum of arr[0...n] = 1
	glp_prob *lp;
	lp = glp_create_prob();
	glp_set_prob_name(lp, "coefficients");
	glp_set_obj_dir(lp, GLP_MIN);
	glp_add_rows(lp, 1);
	glp_set_row_name(lp, 1, "p");
	glp_set_row_bnds(lp, 1, GLP_FX, 1,1);
	glp_add_cols(lp, len);
	int x;
	double ar[50];
	int xar[50];
	int yar[50];
	for (x = 0; x < len; x++) {
		glp_set_col_bnds(lp, x+1, GLP_DB, 0, 1);
		glp_set_obj_coef(lp, x+1, matr[x]);
		ar[x] = 1;
		xar[x] = x+1;
		yar[x] = 1;
	}
	printf("%d\n", len);
	glp_load_matrix(lp, len-1, yar, xar, ar);
	glp_simplex(lp, NULL);
	printf("%f\n", glp_get_obj_val(lp));
	for (x = 0; x < len; x++) {
		arr[x] = glp_get_col_prim(lp, x+1);
	}
	glp_delete_prob(lp);
	printf("=========FINISH======\n");*/
	
	float trial_coefficients[50];
	float current_lowest = 1000;
	//float lowest_coefficients[50];
	
	int c, x;
	float r = 0;
	float cumulative = 0;
	float current_sum = 0;

	/* It's legit I swear */
	/* You idiot, it's the sum of squares which = 1 */
	for (c = 0; c < pow(100, len); c++) {
		cumulative = 0;
		current_sum = 0;
		for (x = 0; x < len-1; x++) {
			r = 1000 - (cumulative * 1000);
			
			trial_coefficients[x] = (((float)(rand()%(2*((int) r)))) / 1000)-1;
			cumulative = cumulative + pow(trial_coefficients[x], 2);
			current_sum += trial_coefficients[x] * matr[x];
		}
		trial_coefficients[len-1] = sqrt(1 - cumulative);
		current_sum += trial_coefficients[len-1] * matr[len-1];
		if (abs(current_sum) < current_lowest) {
			current_lowest = abs(current_sum);
			for (x = 0; x < len; x++) {
				arr[x] = trial_coefficients[x];
			}
		}
		
	}
	//printf("Current Lowest: %f\n", current_lowest);
	int i;
	
}

int main(void) {
	//srand(time(NULL));
	//printf("%f\n", ((float)(rand()%1000))/1000);
	//exit(0);
	int n = 50;

	
	atom atoms[5];
	int num_atoms = 2;
	atoms[0].orbitals[0][0] = 1;
	atoms[0].orbitals[0][1] = 0;
	atoms[0].orbitals[0][2] = 0;
	atoms[0].num_orbitals = 1;
	atoms[0].pos[0] = 0;
	atoms[0].pos[1] = 0;
	atoms[0].pos[2] = 0;
	atoms[0].zeff = 1;
	atoms[1].orbitals[0][0] = 1;
	atoms[1].orbitals[0][1] = 0;
	atoms[1].orbitals[0][2] = 0;
	atoms[1].num_orbitals = 1;
	atoms[1].pos[0] = 1.4;
	atoms[1].pos[1] = 0;
	atoms[1].pos[2] = 0;
	atoms[1].zeff = 1;

	
	/* -1.5 for 1s, -0.375 for 2s, -0.167 for 3s */
	
	float Ss[50][50];
	float Hs[50][50];
	pair xs;
	int x, y, dim = num_atoms, x_orb, y_orb;
	float bond_length = 0.5;

	float distance;
	int a = 0,b = 0;
	for (x = 0; x < num_atoms; x++) {
		for (x_orb = 0; x_orb < atoms[x].num_orbitals; x_orb++) {
			for (y = 0; y < num_atoms; y++) {
				for (y_orb = 0; y_orb < atoms[y].num_orbitals; y_orb++) {
			/// pair calculate_energy(int nA, int lA, int qA, int nB, int lB, int qB, int Za, int Zb, double bond_length, int siz) {
					distance = sqrt(pow(atoms[x].pos[0] - atoms[y].pos[0], 2) + pow(atoms[x].pos[1] - atoms[y].pos[1], 2) + pow(atoms[x].pos[2] - atoms[y].pos[2], 2));
					//printf("DISTANCE %f\n", distance);
					//printf("%d%d%d %f // %d%d%d %f :: %f\n", atoms[x].orbitals[x_orb][0], atoms[x].orbitals[x_orb][1], atoms[x].orbitals[x_orb][2],atoms[x].zeff, atoms[y].orbitals[y_orb][0], atoms[y].orbitals[y_orb][1], atoms[y].orbitals[y_orb][2], atoms[y].zeff, distance * bond_length);
					xs = calculate_energy(atoms[x].orbitals[x_orb][0], atoms[x].orbitals[x_orb][1], atoms[x].orbitals[x_orb][2], atoms[y].orbitals[y_orb][0], atoms[y].orbitals[y_orb][1], atoms[y].orbitals[y_orb][2], atoms[x].zeff, atoms[y].zeff, distance * bond_length, n);
					Ss[a][b] = xs.b;
					if (distance > 0)
						Hs[a][b] = (xs.a/3) + ((atoms[x].zeff * atoms[y].zeff) / (distance * bond_length));
					else
						Hs[a][b] = (xs.a/3);
					//printf("Ss[%d][%d] = %f;; Hs[%d][%d] = %f\n", a, b, Ss[a][b], a, b, Hs[a][b]);
					b++;
				}
			}
			a++;	
			b = 0;		
		}
	}

	printf("H\n");
	for (x = 0; x < a; x++) {
		for (y = 0; y < a; y++) {
			printf("%0.4f\t", Hs[x][y]/1000);
		}
		printf("\n");
	}
	printf("S\n");
	for (x = 0; x < a; x++) {
		for (y = 0; y < a; y++) {
			printf("%0.4f\t", Ss[x][y]/1000);
		}
		printf("\n");
	}
	
	/*if (a != b) {
		printf("ded\n");
		exit(0);
	}*/
	
	float E;
	float prev;
	float curr;
	float m[50][50];
	float switching_points[50];
	int curr_s = 0;
	for (E = -6; E <= 0; E+=0.001) {
		for (x = 0; x < a; x++) {
			for (y = 0; y < a; y++) {
				m[x][y] = Hs[x][y] - E * Ss[x][y];
			}
		}	
		prev = curr;
		curr = det(m, dim);
		if (E != 0) {
			if (switched(prev, curr) == 1) {
				//printf("%f ===================\n", E);
				switching_points[curr_s] = E;
				curr_s++;
				
				
			}
		}
		//printf("%f, %f\n", E, curr);
	}
	float coeff[50];
	float red[50];
	for (x = 0; x < curr_s; x ++) {
		printf("%f - [", switching_points[x]);
		/* Now calculate the coefficients */
		for (y = 0; y < a; y++) {
			//if (y == 0)
				red[y] = Hs[0][y] - (switching_points[x]) * Ss[0][y];
			//else
		//		red[y] = Hs[0][y];
			//printf("%f %f\n", Hs[0][y], Ss[0][y]);
			//printf("Red[%d] %f\n", y, red[y]);
		}
		calculate_coefficients(red, a, coeff);
	//	printf("[");
		for (y = 0; y < a; y++) {
			printf("%f", coeff[y]);
			if (y != a-1)
				printf("\t");
		}
		printf("]\n");
		
	}

	return 1;
}
