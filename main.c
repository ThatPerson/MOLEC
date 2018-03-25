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
	double x, y, z;
} position;

typedef struct {
	int orbitals[50][3];
	int num_orbitals;
	position pos;
	float zeff;
    char symbol;
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

double wavefunction(int n, int l, int q, position read_pos, position pos, float Z) {
	// args is [n, l, m, Z]
	double x = pos.x - read_pos.x;
	double y = pos.y - read_pos.y;
	double z = pos.z - read_pos.z;
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
		printf("%d %d %d\n", n, l, q);
		printf("I don't recognise that orbital..\n");
		exit(1);
	}
	return -1;
}

position arr_pos(double x, double y, double z) {
	position ls;
	ls.x = x;
	ls.y = y;
	ls.z = z;
	return ls;
}

void print_pos(position p) {
    printf("(%f, %f, %f)\n", p.x, p.y, p.z);
    return;
}

double hamiltonian(int nA, int lA, int qA, int nB, int lB, int qB, position posA, position read_pos, float Z, float ZB, position posB) {
/* This is the problem!!! */


	// Numerical Form of hamiltonian
	double alpha = 0.01;
	double beta = 0.001; // Beta needs to be smaller than alpha
	// r component
	//posA[0] = posA[0]/10;
	//posA[1] = posA[1]/10;
	//posA[2] = posA[2]/10;
	double r = sqrt(pow(posA.x, 2) + pow(posA.y, 2) + pow(posA.z, 2));
	double rB = sqrt(pow(posB.x - posA.x, 2) + pow(posB.y - posA.y, 2) + pow(posB.z - posA.z, 2));
	
	double distRA = sqrt(pow(posA.x - read_pos.x, 2) + pow(posA.y - read_pos.y, 2) + pow(posA.z - read_pos.z, 2));
	double distRB = sqrt(pow(posB.x - read_pos.x, 2) + pow(posB.y - read_pos.y, 2) + pow(posB.z - read_pos.z, 2));
	
	
	position dxsquu = arr_pos(posA.x + alpha + beta, posA.y, posA.z);
	position dxsqud = arr_pos(posA.x + alpha - beta, posA.y, posA.z);
	position dxsqdu = arr_pos(posA.x - alpha + beta, posA.y, posA.z);
	position dxsqdd = arr_pos(posA.x - alpha - beta, posA.y, posA.z);
	
	position dysquu = arr_pos(posA.x, posA.y + alpha + beta, posA.z);
	position dysqud = arr_pos(posA.x, posA.y + alpha - beta, posA.z);
	position dysqdu = arr_pos(posA.x, posA.y - alpha + beta, posA.z);
	position dysqdd = arr_pos(posA.x, posA.y - alpha - beta, posA.z);
	
	position dzsquu = arr_pos(posA.x, posA.y, posA.z + alpha + beta);
	position dzsqud = arr_pos(posA.x, posA.y, posA.z + alpha - beta);
	position dzsqdu = arr_pos(posA.x, posA.y, posA.z - alpha + beta);
	position dzsqdd = arr_pos(posA.x, posA.y, posA.z - alpha - beta);
	
	//		dxsquu[0], dxsquu[1], dxsquu[2], dxsquu[3], 
	//		dxsqud[0], dxsqud[1], dxsqud[2], dxsqud[3], 
	//		dxsqdu[0], dxsqdu[1], dxsqdu[2], dxsqdu[3], 
	//		dxsqdd[0], dxsqdd[1], dxsqdd[2], dxsqdd[3]);
	
	// -(2.64589 Power[10, -10]) (D[D[Hx[x, y, z], x], x] - 
   //  D[D[Hx[x, y, z], y], y] - D[D[Hx[x, y, z], z], z]) - 
   //  10/Sqrt[x^2 + y^2 + z^2] Hx[x, y, z];
	
	
	
	/* Work out how to map this onto the actual values.!!!!!!!!!!!!!! */
	
	
	double dxsq = (((wavefunction(nA, lA, qA, read_pos, dxsquu, Z) - wavefunction(nA, lA, qA, read_pos, dxsqud, Z)) / (2*beta)) - ((wavefunction(nA, lA, qA, read_pos, dxsqdu, Z) - wavefunction(nA, lA, qA, read_pos, dxsqdd, Z)) / (2*beta))) / (2*alpha);
	double dysq = (((wavefunction(nA, lA, qA, read_pos, dysquu, Z) - wavefunction(nA, lA, qA, read_pos, dysqud, Z)) / (2*beta)) - ((wavefunction(nA, lA, qA, read_pos, dysqdu, Z) - wavefunction(nA, lA, qA, read_pos, dysqdd, Z)) / (2*beta))) / (2*alpha);
	double dzsq = (((wavefunction(nA, lA, qA, read_pos, dzsquu, Z) - wavefunction(nA, lA, qA, read_pos, dzsqud, Z)) / (2*beta)) - ((wavefunction(nA, lA, qA, read_pos, dzsqdu, Z) - wavefunction(nA, lA, qA, read_pos, dzsqdd, Z)) / (2*beta))) / (2*alpha);
	
	//printf("%lf %lf %lf :SQUARES: %lf %lf %lf\n", posA[0], posA[1], posA[2], dxsq, dysq, dzsq);
	
	double a = 1/2.;
	double b = 1;
	
	//double a = 1;
    double otheratom = 0;
    double potentialenergy = 0;
	//double b = 1;
	double kineticenergy = - a * (dxsq + dysq + dzsq);
    if (distRA != 0)
        potentialenergy = - b * wavefunction(nA, lA, qA, read_pos, posA, Z) * Z / distRA;
    if (rB != 0 && distRB != 0) // Don't count the same atom twice.
        otheratom = - b * wavefunction(nB, lB, qB, read_pos, posB, ZB) * ZB / distRB;

	//printf("%f %f %f\n", kineticenergy, potentialenergy, otheratom);
	double totalenergy = kineticenergy + (potentialenergy+ otheratom);
	
	
	
	//printf("%lf %lf %lf %lf\n", rse, thetase, nA, lA, qAse, chase);
	
	return (double)totalenergy/* * (2.1673 * pow(10, -13))*/;
}

pair calculate_energy(int nA, int lA, int qA, int nB, int lB, int qB, float Za, float Zb, position sA, position sB, int siz) {

	double numerator = 0;
	double denominator = 0;

	double hammy = 0;
	double AWav = 0;
	double BWav = 0;

	position pos;

	double priv_num, priv_den;

	priv_num = 0;
	priv_den = 0;
	double r, theta, phi;
	
	pair xl;
	xl.a = numerator;
	xl.b = denominator;
    
    /*pair calculate_energy(int nA, int lA, int qA, int nB, int lB, int qB, float Za, float Zb, position sA, position sB, int siz) {*/
    position read = arr_pos((sA.x + sB.x) / 2, (sA.y + sB.y) / 2, (sA.z + sB.z) / 2);
    if ((read.x == sA.x && read.y == sA.y && read.z == sA.z) || (read.x == sB.x && read.y == sB.y && read.z == sB.z))
        read = arr_pos(0, 1, 0);
    //read = arr_pos(3, 0, 0);
    float ham = hamiltonian(nA, lA, qA, nB, lB, qB, sA, read, Za, Zb, sB);
    float wav = wavefunction(nA, lA, qA, read, sA, Za);
    xl.a = ham / wav; // Eigenfunction eigenvalue pair.
    
    position A = arr_pos(0, 0, 0);
    position B = arr_pos(2, 0, 0);
    position reads = arr_pos(3, 0, 0);

    float x, y, z, volume_element, step = 0.5, total_sum_oi = 0, oi, total_sum_ham = 0;
    for (x = -30; x <= 100; x+=step) {
        for (y = -10; y <= 10; y += step) {
            for (z = -10; z <= 10; z+=step) {
                volume_element = pow(step, 3);
                oi = wavefunction(nA, lA, qA, arr_pos(x, y, z), sA, Za) * wavefunction(nB, lB, qB, arr_pos(x, y, z), sB, Zb) * volume_element;
                total_sum_oi += oi;
                ham = hamiltonian(nA, lA, qA, nB, lB, qB, sA, arr_pos(x, y, z), Za, Zb, sB) * wavefunction(nB, lB, qB, arr_pos(x, y, z), sB, Zb);
                total_sum_ham += volume_element * ham;
            }
        }
    }
    xl.b = total_sum_oi;
   // xl.a = total_sum_ham;
    xl.a = hamiltonian(nA, lA, qA, nB, lB, qB, sA, read, Za, Zb, sB) / wavefunction(nA, lA, qA, read, sA, Za);
    
    if (sA.x != sB.x) {
        printf("\n\n\nSETTINGS: %d %d %d %d %d %d/// %f %f\n", nA, lA, qA, nB, lB, qB, Za, Zb);
        print_pos(read);
        print_pos(sA);
        print_pos(sB);
        printf("ENERGY: %f \n",hamiltonian(nA, lA, qA, nB, lB, qB, sA, read, Za, Zb, sB) / wavefunction(nA, lA, qA, read, sA, Za));
        printf("XLA %f\n\n\n\n\n\n", xl.a);
    }
    // hamiltonian(1, 0, 0, 1, 0, 0, arr_pos(0, 0, 0), arr_pos(3, 0, 0), 1, 1, arr_pos(6, 0, 0)) / wavefunction(1, 0, 0, arr_pos(3, 0, 0), arr_pos(0, 0, 0), 1)
    
    printf("%f %f\n", total_sum_oi, total_sum_ham);
    
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

int read_settings(atom *atoms, char * filename) {

	atoms[0].pos.x = 0;
	atoms[0].pos.y = 0;
	atoms[0].pos.z = 0;

	//atoms[1].pos[0] = 1.4;
	atoms[1].pos.x = 7;
	atoms[1].pos.y = 0;
	atoms[1].pos.z = 0;
    
    

    FILE* fp;
    char buffer[255];

    fp = fopen(filename, "r");
    int mode;
    int orbital;
    int current_atom = 0;
    int buff = 0;
    int cu = 0;
    char tmp_str[50];
    float position[3];
    int curr_p = 0;
    while(fgets(buffer, 255, (FILE*) fp)) {
        int i;
        mode = 0;
        orbital = 0;
        cu = 0;
        strcpy(tmp_str, "");
        buff = 0;
        curr_p = 0;
        for (i = 0; i < 255; i++) {
            if (buffer[i] == 0 || buffer[i] == 10) 
                break;
            switch (mode) {
                case 0: mode++; atoms[current_atom].symbol = buffer[i]; break;
                case 1: if (buffer[i] == 91) mode++; break;
                case 2: 
                    if (buffer[i] == 93) {
                        mode++; 
                        tmp_str[buff] = '\0';
                        atoms[current_atom].zeff = atof(tmp_str);
                    } else if (buffer[i] >= '0' && buffer[i] <= '9') {
                        tmp_str[buff] = buffer[i];
                        buff++;
                    }
                    break;
                case 3: 
                    if (buffer[i] == 40) {
                        mode++; 
                        atoms[current_atom].num_orbitals = 0;
                        cu = 0;
                    } 
                    break;
                case 4: 
                    if (buffer[i] == 44) {
                        atoms[current_atom].num_orbitals++;
                        cu=0;
                    } 
                    if (buffer[i] >= '0' && buffer[i] <= '9') {
                        atoms[current_atom].orbitals[atoms[current_atom].num_orbitals][cu] = buffer[i] - '0';
                        cu++;
                    }
                    if (buffer[i] == 41) { 
                        mode++; 
                        atoms[current_atom].num_orbitals++; 
                    }
                    break;
                case 5:
                    if (buffer[i] == 123) {
                        curr_p = 0;
                        mode++;
                        strcpy(tmp_str, "");
                        buff = 0;
                    }
                    break;
                case 6:
                    if (buffer[i] == 44) {
                        tmp_str[buff] = '\0';
                        position[curr_p] = atof(tmp_str);
                        strcpy(tmp_str, "");
                        buff = 0;
                        curr_p++;
                        if (curr_p > 2) {
                            printf("Out of bounds.\n");
                            break;
                        }
                    
                    } else if (buffer[i] >= '0' && buffer[i] <= '9') {
                        tmp_str[buff] = buffer[i];
                        buff++;
                    } else if (buffer[i] == 125) {
                        atoms[current_atom].pos = arr_pos(position[0], position[1], position[2]);
                        mode++;
                    }
                    break;
                    
                    
                    
                    
                    
            } 
        }
        printf("THIS ATOM;;; \n%c -- %f   %d\n", atoms[current_atom].symbol, atoms[current_atom].zeff, atoms[current_atom].num_orbitals);

        for (i = 0; i < atoms[current_atom].num_orbitals; i++) {
            printf("\t%d %d %d\n", atoms[current_atom].orbitals[i][0], atoms[current_atom].orbitals[i][1], atoms[current_atom].orbitals[i][2]);
        }
        printf("at position (%f, %f, %f)\n", atoms[current_atom].pos.x, atoms[current_atom].pos.y, atoms[current_atom].pos.z);
        current_atom++;
    }

    fclose(fp);


    
    return current_atom;
}

int main(int argc, char*argv[]) {
	//srand(time(NULL));
	//printf("%f\n", ((float)(rand()%1000))/1000);
	//exit(0);
	int n = 50;

	
	atom atoms[5];
    if (argc < 2) {
        printf("Please pass a filename\n");
        return 1;
    }
	int num_atoms = read_settings(atoms, argv[1]);
	

    

	double bonding[100];
	double antibonding[100];
	double nonbonding[100];
	double psd = 4;
    int current = 0;
    /*for (psd = 0.5; psd < 50; psd += 0.5) { // If the upper limit is over half that of the total range scanned it fails?

        printf("LSLSLS %f\n", psd);
        atoms[0].pos.x = 0;
        atoms[1].pos.x = psd;*/
	
        float Ss[50][50];
        float Hs[50][50];
        pair xs;
        int x, y, dim = num_atoms, x_orb, y_orb;
        float bond_length = 0.5;

        float distance = 0;
        int a = 0,b = 0;
        for (x = 0; x < num_atoms; x++) {
            for (x_orb = 0; x_orb < atoms[x].num_orbitals; x_orb++) {
                for (y = 0; y < num_atoms; y++) {
                    
                    for (y_orb = 0; y_orb < atoms[y].num_orbitals; y_orb++) {
                        xs = calculate_energy(atoms[x].orbitals[x_orb][0], atoms[x].orbitals[x_orb][1], atoms[x].orbitals[x_orb][2], atoms[y].orbitals[y_orb][0], atoms[y].orbitals[y_orb][1], atoms[y].orbitals[y_orb][2], atoms[x].zeff, atoms[y].zeff, atoms[x].pos, atoms[y].pos, n);
                        Ss[a][b] = xs.b;
                        Hs[a][b] = (xs.a);
                        Hs[b][a] = (xs.a);
                        Ss[b][a] = (xs.b);

                        b++;
                    }
                }
                a++;	
                b = 0;		
            }
        }

        printf("\tH\n\t\t");
        for (x = 0; x < a; x++) {
            for (y = 0; y < a; y++) {
                printf("%0.4f\t", Hs[x][y]);
            }
            printf("\n\t\t");
        }
        printf("\n\tS\n\t\t");
        for (x = 0; x < a; x++) {
            for (y = 0; y < a; y++) {
                printf("%0.4f\t", Ss[x][y]);
            }
            printf("\n\t\t");
        }
        printf("\n");
        float E;
        float prev;
        float curr;
        float m[50][50];
        float switching_points[50];
        int curr_s = 0;
        for (E = -6; E <= 6; E+=0.001) {
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
        bonding[current] = switching_points[0];
        antibonding[current] = switching_points[1];
        current++;
   /* }
    int i;
	for (i = 0; i < current; i++) {
		printf("%f, %f, %f\n", (0.5 * i) + 0.5, bonding[i], antibonding[i]);
	}*/
	return 1;
}
