#include <stdio.h>
#include <math.h>

// Most commonly adjusted paramaters
#define dimensions 3				// maximum 4, assuming a span of 64
#define span 8						// width of lattice along every dimension
#define Tmin 1.0					// minimum value of beta
#define Tmax 10.0					// maximum value of beta
#define samples 200					// # of samples between betamin and betamax
#define	experiments 1				// # of separate experiments to compile
#define runcount 100000				// # of averaging runs at each value of beta
#define runsteps 1000				// # of Monte Carlo steps between each averaging run

// Other paramaters and constants
#define coldstart 0					// -1 = coldstart, 0 = hotstart
#define equlibriate 100000			// # of Monte Carlo steps given to reach equlibrium at each value of beta
#define Boltzmann 1.0				// Boltzmann constant (used only for scaling heat capacity)

// Global variables
double beta;						// Thermodynamic beta  = 1 / T
long sample_critical;				// sample number of critical beta
long modulus;						// number of cells in the lattice
char lattice[16777216];				// allocation of memory for the lattice, 1 byte per cell
									//	sufficient space for 64 * 64 * 64 * 64 lattice size
long neighbours[2 * dimensions];	// list of nearest neighbours for a given cell
long reach[dimensions];				// required reach to nearest neighbours
double energy, energy2, magnetiz;	// statistics, averages per cell
double energy_list[samples];		// results arrays holding the output at each sample point for a single experiment
double energy2_list[samples];
double magnetiz_list[samples];
double Temp_list[samples];
double heatcapacity_list[samples];
long critical_histogram[samples];	// histogram of the location of critical beta (sample no.) over many experiments

// Function declarations
unsigned long rnd(void);			// Return a random integer x between 0 and 1<<31
double rnd_double(void);			// Return a random double 0 <= x < 1
unsigned long rnd_long(long);		// Return a random integer 0 <= x < n
double abs_double(double);			// abs for double precision numbers
void newlattice(void);				// Prepare a new lattice
void render(void);					// Display the lattice
long rnd_cell(void);				// Pick a random cell in the lattice
void find_neighbours(long n);		// Find the nearest neighbours of the cell
long sum_neighbours(void);			// Sum the states of the neighbours
long cell_energy(long n);			// Sum the interaction energy attributable to cell n
void flip(long);					// Flip cell n
void stats(void);					// Calculate the stats for the lattice 
void experiment(void);				// Run a single experiment from betamin to betamax
void metropolis(void);				// Select and flip a cell according to Metropolis
void print_stats(void);				// Output the latest stats
void fprint_experiment(void);		// Output the results of a simulation run
void fprint_histogram(void);		// Output the histogram of critical values of Beta

int main()
{
	long i;

	for (i = 0; i < samples; i++)
	{
		critical_histogram[i] = 0;
	}

	for (i = 0; i < experiments; i++)	
	{
		printf("%d...",i);
		experiment();
		critical_histogram[sample_critical] += 1;
	}

	fprint_experiment();
	fprint_histogram();

	printf("\nSimulation complete. Press any key to exit.\n");

	getchar();
	return 0;
}

void experiment(void)
{
	long i, j, k;
	double Temp, Tstep, b, h, e, e2;
	double h_max = 0.0;

	newlattice();								
	Tstep = (Tmax - Tmin) / (samples - 1);
	if (!coldstart)								
	{
		Temp = Tmax;
		Tstep = - Tstep;
	}
	else
		Temp = Tmin;

	// zero the averaging arrays
	for (i = 0; i < samples; i++)				
	{
		energy_list[i] = 0.0;
		energy2_list[i] = 0.0;
		magnetiz_list[i] = 0.0;
	}

	// iterate over samples, runs, and Monte-Carlo steps
	for (i = 0; i < samples; i++)
	{
		//printf("\tTemp = %f\n",Temp);
		Temp_list[i] = Temp; beta = 1.0 / Temp;

		for (k = 0; k < equlibriate; k++)
				metropolis();
		for (j = 0; j < runcount; j++)
		{
			for (k = 0; k < runsteps; k++)
				metropolis();
			stats();
			energy_list[i] += energy;
			energy2_list[i] += energy2;
			magnetiz_list[i] += abs_double(magnetiz);
		}
		Temp += Tstep;
	}

	// complete averaging and calculate heat capacities
	for (i = 0; i < samples; i++)
	{
		b = 1.0 / Temp_list[i];
		e = energy_list[i] / runcount;
		e2 = energy2_list[i] / runcount;

		magnetiz_list[i] /= runcount;
		energy_list[i] = e;
		energy2_list[i] = e2;
		h = Boltzmann * b *  b * (e2 - e * e);
		heatcapacity_list[i] = h;

		// screening for the critical point
		if (h > h_max)
		{
			h_max = h;
			sample_critical = i;
		}
	}

}
		
// Select and flip a cell according to Metropolis
void metropolis(void)				
{
	long n, e;
	double p;

	n = rnd_cell();
	e = cell_energy(n);
	if (e < 0)					// a flip will raise the energy of the system
	{
		p = exp(beta * e * 2.0);
		if (rnd_double() < p)
		{
			flip(n);
		}
	}
	else
		flip(n);

}

// Data export routines

// Output the latest stats
void print_stats(void)			
{
	printf(" <magnetization> %f\n", magnetiz);
	printf(" <energy> %f\n", energy);
	printf(" <energy^2> %f\n", energy2);
}

// Output the results of a simulation run
void fprint_experiment(void)			
{
	FILE *fp;
	int i;

	fopen_s(&fp, "experiment.txt","w");
	if (fp != NULL)
	{
		fprintf(fp,"Temp\tMagnetization\tHeat capacity\tEnergy\n");
		for (i = 0; i < samples; i++)
			fprintf(fp,"%f\t%f\t%E\t%f\n", Temp_list[i], magnetiz_list[i], heatcapacity_list[i], energy_list[i]);
	fclose(fp);	
	}

	fopen_s(&fp, "signature.txt","w");
	if (fp != NULL)
	{
		fprintf(fp,"Dimensions = %d\n",dimensions);
		fprintf(fp,"Span = %d\n",span);
		fprintf(fp,"No. steps to equlibriate at each temp = %d\n",equlibriate);
		fprintf(fp,"No. averages at each temp = %d\n",runcount);
		fprintf(fp,"No. step between each average = %d\n",runsteps);
		fprintf(fp,"No. experiments = %d\n",experiments);
		fclose(fp);	
	}
}

// Output the histogram of critical values of Beta
void fprint_histogram(void)				
{
	FILE *fp;
	int i;

	fopen_s(&fp, "criticaltemp.txt","w");
	if (fp != NULL)
	{
		fprintf(fp,"Sample\tTemp\tCount\n");
		for (i = 0; i < samples; i++)
			fprintf(fp,"%d\t%f\t%d\n", i, Temp_list[i], critical_histogram[i]);
	fclose(fp);	
	}
}

// Utility routines
unsigned long rnd(void)				// Return a random integer x between 0 and 1<<31
{
	static unsigned long X;
	return X = (unsigned long) ((unsigned long long)X * 3141592621ULL + 1ULL);
	// use the 64 integer type ULL to avoid premature truncation
	// the modulus is provided by the typecast (unsigned long)
}

double rnd_double(void)
{
	return ((double) rnd() / 4294967295.0);
}


unsigned long rnd_long(long n)		// Return a random integer 0 <= x < n
{
	return (unsigned long) (rnd_double() * (double) n);
}

double abs_double(double x)			// abs for double precision numbers
{
	if (x > 0)
		return x;
	else
		return -x;
}

// Ising model routines
void newlattice(void)				// Prepare a new lattice
{
	long i;

	// calculate the total size of the lattice
	modulus = (long) pow((double) span, dimensions);

	// set up the array of reach, indicating the offsets to the nearest neighbour in each dimension
	for (i = 0; i<dimensions; i++)
		reach[i] = (long) pow((double) span, i);
	
	// prepare the lattice with either hot or cold start
	for (i=0; i<modulus; ++i)
	{
		if (coldstart || (rnd() & (1<<31)))
			lattice[i] = 2;
		else
			lattice[i] = 0;
	}
}

long rnd_cell(void)					// Pick a random cell in the lattice
{
	return rnd_long(modulus);
}

void find_neighbours(long n)		// Find the nearest neighbours of the given cell
{
	int i;

	for (i = 0; i<dimensions; i++)
	{
		// in each dimension, i, there are two neighbours at a distance reach[i] away
		neighbours[2*i] = (n + reach[i] + modulus) % modulus;		// make sure to bring the offset back into range
		neighbours[2*i+1] = (n - reach[i] + modulus) % modulus;
	}
}

long sum_neighbours(void)			// Sum the states of the neighbours
{
	long i;
	long e=0;

	for (i=0; i < 2*dimensions; i++)
		e += lattice[neighbours[i]] - 1;

	return e;
}

long cell_energy(long n)			// Sum the interaction energy attributable to cell n
{
	find_neighbours(n);
	return -1 *(lattice[n] - 1) * sum_neighbours();
}

void flip(long n)					// Flip cell n
{
	lattice[n] = 2 - lattice[n];
}

void stats(void)					// Calculate the stats for the lattice 
{
	long i;
	double sum_e = 0.0;
	double sum_m = 0.0;
	double z = modulus * 2.0;		// normalization factor includes 2.0 since each interaction energy has been counted twice
	
	for(i=0; i<modulus; i++)
	{
		sum_m += lattice[i] - 1;
		sum_e += cell_energy(i);
	}

	energy = sum_e / z;						
	energy2 = (sum_e *sum_e) / (z * z);
	magnetiz = sum_m / modulus;

}

void render(void)					// Display the lattice
{
	long i, j, n;
	n = 0;

	for (j=0; j<span; ++j)
	{
		for (i=0; i<span; ++i)
			if (lattice[n++] == 2)
				printf("+ ");
			else
				printf("- ");
		printf("\n");
	}
}

// Helpful resource
// http://msdn.microsoft.com/en-us/library/fw5abdx6.aspx
