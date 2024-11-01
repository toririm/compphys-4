#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define LOG_FILENAME "LJ.csv"	// output file for log data
#define XYZ_FILENAME "LJ.xyz"	// output file for structure data
#define LOG_STEP 10 			// output step intercal of log data
#define XYZ_STEP 10				// output step interval of xyz data
#define PRINT_STEP 1000			// output step interval of time
#define N 256					// total number of atoms
#define ROW 4					// 4*ROW^3 = N
#define BOXLX 100.0				// x-side length of unit cell [angstrome]
#define BOXLY 100.0				// y-side length of unit cell [angstrome]
#define BOXLZ 100.0				// z-side length of unit cell [angstrome]
#define SIGMA 3.405				// Lennard-Jones parameter [angstrome]
#define EPSILON 0.01032			// Lennard-Jones parameter [eV]
#define CUTOFF 3.0				// cutoff length [SIGMA]
#define MASS 39.95				// mass of Ar [amu]
#define T_INIT 300.0			// initial temperature [K]
#define TIME_END 100.0			// time length of simulation [ps]
#define TIME_STEP 0.01			// time step [ps]
#define RAND_SEED 10			//seed of random number generator
#define ps 1.0e-12				// 1 pico second [s]
#define angstrom 1.0e-10		// 1 angstrome [m]
#define eV 1.602176e-19			// 1 electron volt [J]
#define amu 1.660539e-27		// 1 atomic mass unit [kg]
#define kb 8.161734e-5				// Boltzmann constant [eV/K]

void calc_LJ_potential_force(double *fx, double *fy, double *fz,
							 double *rx, double *ry, double *rz,
							 double *PE);

int main(int argc, char* argv[])
{
	printf("===============================\n");
	printf("=  M D   S I M U L A T I O N  =\n");
	printf("=  LJ-NVE  v1.0 (2024.10.16)  =\n");
	printf("===============================\n");
	clock_t start_t = clock();
	double rx[N],  ry[N],  rz[N];		// coordinates
	double r0x[N], r0y[N], r0z[N];		// coordinates at the first step
	double vx[N],  vy[N],  vz[N];		// velocity
	double fx[N],  fy[N],  fz[N];		// force
	double f0x[N], f0y[N], f0z[N];		// force at the previous step
	double m;							// mass of atoms
	double rij;							// inter-atomic distance
	double dxij, dyij, dzij;			// relative position vector
	double msdx, msdy, msdz, msd;		// mean-square displacement
	double t, dt;						// time, time step
	double KE;							// kinetic energy [eV]
	double PE;							// potential energy [eV]
	double latticeconst;				// lattice constant at initial state
	double massunitconv;				// unit conversion factor for mass
	double T;							// temperature [K]
	int counter;						// MD simulation step counter
	int i,j,k,l;
	FILE *fp1, *fp2, *fp3;

	fp1 = fopen(LOG_FILENAME, "w");
	if(fp1 == NULL){
		printf("ERROR: File open failed.\n"); exit(EXIT_FAILURE);
	}
	fp2 = fopen(XYZ_FILENAME, "w");
	if(fp2 == NULL){
		printf("ERROR: File open failed.\n"); exit(EXIT_FAILURE);
	}
	fprintf(fp1,"t,KE,PE,KE+PE,T,MSD\n");
	if(CUTOFF*SIGMA > 0.5*BOXLX || CUTOFF*SIGMA > 0.5*BOXLY || CUTOFF*SIGMA > 0.5*BOXLZ){
		printf("ERROR: Unitcell size is too small.\n");
		exit(EXIT_FAILURE);
	}

	// Initialize
	dt = TIME_STEP;
	massunitconv = amu * angstrom * angstrom / ps / ps / eV;
	m = MASS * massunitconv;
	srand(RAND_SEED);
	// Initial configuration (face-centered cubic)
	latticeconst = pow(2.0, 1.0/6.0) * SIGMA * sqrt(2.0);
	printf("\n### PARAMETERS ###\n");
	printf("mass:\t\t\t%12.6f [amu]\n", MASS);
	printf("T_init:\t\t\t%12.6f [K]\n", T_INIT);
	printf("timestep:\t\t%12.6f [ps]\n", dt);
	printf("simulation time:\t%12.6f [ps]\n", TIME_END);
	printf("LJ_sigma:\t\t%12.6f [angstrom]\n", SIGMA);
	printf("LJ_epsilon:\t\t%12.6f [eV]\n", EPSILON);
	printf("##################\n\n");
	i = 0;
	printf("Generating Coordinates (Ar fcc)...\n");
	for(j=0; j<ROW; ++j){
		for(k=0;k<ROW; ++k){
			for(l=0; l<ROW; ++l){
				rx[i] = latticeconst * (j - 0.5 * ROW) + 0.5 * BOXLX;
				ry[i] = latticeconst * (k - 0.5 * ROW) + 0.5 * BOXLY;
				rz[i] = latticeconst * (l - 0.5 * ROW) + 0.5 * BOXLZ;
				r0x[i] = rx[i];
				r0y[i] = ry[i];
				r0z[i] = rz[i];
				++i;
				rx[i] = latticeconst * (j + 0.5 - 0.5 * ROW) + 0.5 * BOXLX;
				ry[i] = latticeconst * (k + 0.5 - 0.5 * ROW) + 0.5 * BOXLY;
				rz[i] = latticeconst * (l - 0.5 * ROW) + 0.5 * BOXLZ;
				r0x[i] = rx[i];
				r0y[i] = ry[i];
				r0z[i] = rz[i];
				++i;
				rx[i] = latticeconst * (j + 0.5 - 0.5 * ROW) + 0.5 * BOXLX;
				ry[i] = latticeconst * (k - 0.5 * ROW) + 0.5 * BOXLY;
				rz[i] = latticeconst * (l + 0.5 - 0.5 * ROW) + 0.5 * BOXLZ;
				r0x[i] = rx[i];
				r0y[i] = ry[i];
				r0z[i] = rz[i];
				++i;
				rx[i] = latticeconst * (j - 0.5 * ROW) + 0.5 * BOXLX;
				ry[i] = latticeconst * (k + 0.5 - 0.5 * ROW) + 0.5 * BOXLY;
				rz[i] = latticeconst * (l + 0.5 - 0.5 * ROW) + 0.5 * BOXLZ;
				r0x[i] = rx[i];
				r0y[i] = ry[i];
				r0z[i] = rz[i];
				++i;
			}
		}
	}
	printf("Generating Velocities (Box-Muller)...\n");
	// set initial velocity by the Box-Muller method
	for(i=0; i<N; ++i){
		vx[i] = sqrt(-2.0*kb*T_INIT/m*log((double)rand()/RAND_MAX))
			    * cos(2.0*M_PI*(double)rand()/RAND_MAX);
		vy[i] = sqrt(-2.0*kb*T_INIT/m*log((double)rand()/RAND_MAX))
			    * cos(2.0*M_PI*(double)rand()/RAND_MAX);
		vz[i] = sqrt(-2.0*kb*T_INIT/m*log((double)rand()/RAND_MAX))
			    * cos(2.0*M_PI*(double)rand()/RAND_MAX);
	}
	printf("Completed!\n\n");
	// initial force calculation
	calc_LJ_potential_force(fx, fy, fz, rx, ry, rz, &PE);

	printf("Start MD simulation\n");
	counter = 0;
	printf("Time=%7.1f ps", 0.0);
	for(t = 0.0; t <= TIME_END; t += dt){
		// Kinetic energy calculation
		KE = 0.0;
		for(i=0; i<N; ++i){
			KE += 0.5 * m * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
		}
		// Temperature calculation
		T = 2.0 / (3.0 * kb) * KE / (double)N;
		
		// MSD calculation
		if(counter % LOG_STEP == 0){
			for(i=0; i<N; ++i){
				msdx = rx[i] - r0x[i];
				msdy = ry[i] - r0y[i];
				msdz = rz[i] - r0z[i];
				if(msdx > 0.5 * BOXLX)	msdx -= BOXLX;
				if(msdx <-0.5 * BOXLX)	msdx += BOXLX;
				if(msdy > 0.5 * BOXLY)	msdy -= BOXLY;
				if(msdy <-0.5 * BOXLY)	msdy += BOXLY;
				if(msdz > 0.5 * BOXLZ)	msdz -= BOXLZ;
				if(msdz <-0.5 * BOXLZ)	msdz += BOXLZ;
				msd += msdx * msdx + msdy * msdy + msdz * msdz;
			}
			msd = msd / N;
		}

		if(counter % LOG_STEP == 0){
			// print kinetic, potential and total energies
			fprintf(fp1, "%f,%f,%f,%f,%f,%f\n", t, KE, PE, KE + PE, T, msd);
		}
		if(counter % XYZ_STEP == 0){
			// print structure
			fprintf(fp2, "%d\n", N);
			fprintf(fp2, "Lennard-Jones potential\n");
			for(i=0; i<N; ++i){
				fprintf(fp2, "Ar %f %f %f\n", rx[i], ry[i], rz[i]);
			}
		}
		if(counter % PRINT_STEP == 0){
			printf("\rTime=%7.1f ps", t);
			fflush(stdout);
		}
		// velocity-Verlet algorithm
		for(i=0; i<N; ++i){
			rx[i] += vx[i] * dt + 0.5 * fx[i] / m * dt * dt;
			ry[i] += vy[i] * dt + 0.5 * fy[i] / m * dt * dt;
			rz[i] += vz[i] * dt + 0.5 * fz[i] / m * dt * dt;
			f0x[i] = fx[i];
			f0y[i] = fy[i];
			f0z[i] = fz[i];
		}
		calc_LJ_potential_force(fx, fy, fz, rx, ry, rz, &PE);
		for(i=0; i<N; ++i){
			vx[i] += 0.5 * (fx[i] + f0x[i]) / m * dt;
			vy[i] += 0.5 * (fy[i] + f0y[i]) / m * dt;
			vz[i] += 0.5 * (fz[i] + f0z[i]) / m * dt;
		}
		// Periodic boundary condition
		for(i=0; i<N; ++i){
			if(rx[i] > BOXLX) 	rx[i] -= BOXLX;
			if(rx[i] < 0.0)		rx[i] += BOXLX;
			if(ry[i] > BOXLY) 	ry[i] -= BOXLY;
			if(ry[i] < 0.0)		ry[i] += BOXLY;
			if(rz[i] > BOXLZ) 	rz[i] -= BOXLZ;
			if(rz[i] < 0.0)		rz[i] += BOXLZ;
		}
		counter++;
	}
	printf("\nCompleted!\n\n");
	fclose(fp1);
	fclose(fp2);
	clock_t end_t = clock();
	double elapsed_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
	printf("Elapsed time: %f [sec]\n", elapsed_t);
	return EXIT_SUCCESS;
}

void calc_LJ_potential_force(double *fx, double *fy, double *fz,
							 double *rx, double *ry, double *rz,
							 double *PE)
{
	int i, j;
	double rij, fij;			// distance, force
	double dxij, dyij, dzij;	// relative position vector

	(*PE) = 0.0;
	for(i=0; i<N; ++i){
		fx[i] = 0.0; fy[i] = 0.0; fz[i] = 0.0;
	}
	for(i=0; i<N; ++i){
		for(j= i+1; j<N; ++j){
			dxij = rx[i] - rx[j];
			dyij = ry[i] - ry[j];
			dzij = rz[i] - rz[j];

			// Minimum image convention
			if(dxij > 0.5 * BOXLX)	dxij -= BOXLX;
			if(dxij <-0.5 * BOXLX)	dxij += BOXLX;
			if(dyij > 0.5 * BOXLY)	dyij -= BOXLY;
			if(dyij <-0.5 * BOXLY)	dyij += BOXLY;
			if(dzij > 0.5 * BOXLZ)	dzij -= BOXLZ;
			if(dzij <-0.5 * BOXLZ)	dzij += BOXLZ;

			rij = sqrt(dxij * dxij + dyij * dyij + dzij * dzij);
			if(rij < CUTOFF * SIGMA){
				(*PE) += 4.0 * EPSILON * (pow(SIGMA/rij,12.0)-pow(SIGMA/rij,6.0));
				fij   = -4.0 * EPSILON * (-12.0 * pow(SIGMA/rij, 12.0) / rij
										  + 6.0 * pow(SIGMA/rij, 6.0) / rij);
				fx[i] += fij * dxij / rij;
				fy[i] += fij * dyij / rij;
				fz[i] += fij * dzij / rij;
				fx[j] -= fij * dxij / rij;
				fy[j] -= fij * dyij / rij;
				fz[j] -= fij * dzij / rij;
			}
		}
	}
	return;
}
