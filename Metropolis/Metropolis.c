#include <stdio.h>
#include <math.h>

// Helpful resource
// http://msdn.microsoft.com/en-us/library/fw5abdx6.aspx

// Preprocessor numeric constants
#define dimensions 2				// maximum 4, assuming a span of 64
#define span 64						// width of lattice along every dimension
#define coldstart 1					// 0 = coldstart, 1 = hotstart
#define betamax	1.0					// maximum value of beta
#define betamin 0.0					// minimum value of beta
#define samples 200					// # of steps in beta across which to sample the phenomena
#define runcount 1000				// # of runs for averaging at each temperature
#define equlibriate 10000			// no. steps to equlibriate
	
// Global variables
double beta;						// Thermodynamic beta  = 1 / T
long modulus;						// number of cells in the lattice
char lattice[16777216];				// allocation of memory for the lattice, 1 byte per cell
									//	sufficient space for 64 * 64 * 64 * 64 lattice size
long neighbours[2 * dimensions];	// list of nearest neighbours for a given cell
long reach[dimensions];				// required reach to nearest neighbours
double energy, energy2, magnetiz;	// statistics, averages per cell
double energy_list[samples];		// results arrays holding the output at each sample point
double energy2_list[samples];
double magnetiz_list[samples];
double beta_list[samples];
double heatcapacity_list[samples];

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
void metropolis(void);				// Select and flip a cell according to Metropolis
void stats(void);					// Calculate the stats for the lattice 
void print_results(void);			// Output the latest stats
void fprint_results(void);			// Output the results of a simulation run

// Library routines
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

// Metropolis algorithm implementation
int main()
{
	long i, j, k;
	double beta_step;

	// initialize
	newlattice();								// initialize a lattice

	beta_step = (betamax - betamin) / (samples - 1);
	if (coldstart)								// beta starting point and step sign
	{
		beta = betamax;
		beta_step = - beta_step;
	}
	else
		beta = betamin;

	for (i = 0; i < samples; i++)				// zero the averaging arrays
	{
		energy_list[i] = 0.0;
		energy2_list[i] = 0.0;
		magnetiz_list[i] = 0.0;
	}

	// iterate over samples, runs, and Monte-Carlo steps
	for (i = 0; i < samples; i++)
	{
		printf("Beta = %f\n",beta);
		beta_list[i] = beta;
		for (j = 0; j < runcount; j++)
		{
			for (k = 0; k < equlibriate; k++)
				metropolis();
			stats();
			energy_list[i] += energy;
			energy2_list[i] += energy2;
			magnetiz_list[i] += abs_double(magnetiz);
		}
		beta += beta_step;
	}

	// complete averaging and calculate heat capacity
	
	for (i = 0; i < samples; i++)
	{
		energy_list[i] /= runcount;
		energy2_list[i] /= runcount;
		magnetiz_list[i] /= runcount;
		heatcapacity_list[i] = beta_list[i] * beta_list[i] * (energy2_list[i] - energy_list[i] * energy_list[i]);
	}

	printf("Simulation complete. Press any key to exit.\n");
	fprint_results();
	getchar();
	return 0;
}

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
			
void metropolis(void)				// Select and flip a cell according to Metropolis
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

void stats(void)					// Calculate the stats for the lattice 
{
	long i;
	double e;
	double sum_e = 0.0;
	double sum_e2 = 0.0;
	double sum_m = 0.0;
	
	for(i=0; i<modulus; i++)
	{
		sum_m += lattice[i] - 1;
		e = cell_energy(i);
		sum_e += e;
		sum_e2 += e*e;
	}

	energy = sum_e / modulus;
	energy2 = sum_e2 / modulus;
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

void print_results(void)			// Output the latest stats
{
	printf(" <magnetization> %f\n", magnetiz);
	printf(" <energy> %f\n", energy);
	printf(" <energy^2> %f\n", energy2);
}

void fprint_results(void)			// Output the results of a simulation run
{
	FILE *fp;
	int i;

	fopen_s(&fp, "metro.txt","w");
	if (fp != NULL)
	{
		fprintf(fp,"beta\t Magnetization\t Heat capacity\t Energy\n");
		for (i = 0; i < samples; i++)
			fprintf(fp,"%f\t%f\t%f\t%f\n",beta_list[i],magnetiz_list[i],heatcapacity_list[i],energy_list[i]);
	fclose(fp);	
	}

}