#include <stdio.h>
#include <math.h>

// Helpful resource
// http://msdn.microsoft.com/en-us/library/fw5abdx6.aspx

// Preprocessor numeric constants
#define dimensions 2				// maximum 4, assuming a span of 64
#define span 16						// width of lattice along every dimension
#define coldstart -1					// -1 = coldstart, 0 = hotstart
#define betamax	0.5					// maximum value of beta
#define betamin 0.3					// minimum value of beta
#define samples 200					// # of samples between betamin and betamax
#define	experiments 1				// # of separate experiments to compile
#define Boltzmann 1.0				// Boltzmann constant (used only for scaling heat capacity)
#define runsteps 100000				// number of Monte Carlo steps between each check of the histogram
#define	referencesteps 1000000		// number of Monte Carlo steps to establish the reference histogram
#define reference_level 1000		// minimum count in the reference histogram to include an energy level
#define flatness_criterion 0.80		// criterion for testing the flatness of the histogram
#define iterationlimit 1000			// criterion for avoiding stuck random walks
#define ln_f_initial 1.0			// initial value of the adjustment factor
#define ln_f_final .0001			// final value of the adjustment factor
#define ln_2 0.6931471806			// ln_g[lowest energy configuration] = ln(2), for normalization

// Macro for addressing the density of states array with an energy level
#define indx(E)	(E + 2 * modulus) / 4
	
// Global variables
double beta;						// Thermodynamic beta  = 1 / T
double beta_critical;				// Beta at the critical point
long sample_critical;				// sample number of critical beta
long modulus;						// number of cells in the lattice
char lattice[16777216];				// allocation of memory for the lattice, 1 byte per cell
double ln_g[16777217];				// allocation of memory for the density of states, as natural logarithm
long hist[16777217];				// allocation of memory for the histogram
long hist_ref[16777217];			// allocation of memory for the reference histogram
double magnetiz_avg[16777217];		// allocation of memory for the average magnetization at each energy level
									//	sufficient space for 64 * 64 * 64 * 64 lattice size
double ln_f;						// adjustment factor for density of states, as natural logarithm
long neighbours[2 * dimensions];	// list of nearest neighbours for a given cell
long reach[dimensions];				// required reach to nearest neighbours
long lattice_energy;				// energy of the current lattice
long lattice_magnetiz;				// magnetization of the current lattice
double energy, energy2, magnetiz;	// statistics, averages per cell
double energy_list[samples];		// results arrays holding the output at each sample point
double energy2_list[samples];
double magnetiz_list[samples];
double beta_list[samples];
double heatcapacity_list[samples];
long critical_histogram[samples];	// histogram of the location of critical beta (sample no.) over many experiments
double	critical_list[experiments];	// list of the critical beta values found at each experiment

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
void print_stats(void);				// Output the latest stats
void fprint_experiment(void);		// Output the results of a simulation run
void zero_dos (void);				// initialize density of states by setting ln_g[i] = 0 (=> g[i] = 1)
void zero_hist (void);				// initialize histogram and the array used for the calculation of average magnetization
void wanglandau(void);				// Select a trial cell and transition the lattice according to Wang Landau
long flat(void);					// check the flatness of the histogram against the criterion
void normalize_dos(void);			// normalize the density states so that ln[ground state] = ln(2)
void num_sums(void);				// perform the numerical sums over the density of states
void microaverage(void);			// calculate the microcanonical average magnetization at each energy level
int experiment(void);				// Run a single experiment from betamin to betamax
void fprint_histogram(void);		// Output the histogram of critical values of Beta
void zero_hist_ref(void);			// initialize reference histogram 

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
	long i;

	for (i = 0; i < samples; i++)				// zero the histogram
	{
		critical_histogram[i] = 0;
	}

	for (i = 0; i < experiments; i++)				// zero the histogram
	{
		printf("%d...",i);
		while (! experiment());
		critical_histogram[sample_critical] += 1;
	}

	fprint_experiment();
	fprint_histogram();

	printf("\nSimulation complete. Press any key to exit.\n");

	getchar();
	return 0;
}

int experiment()
{
	long i;
	long n;

	// initialize
	newlattice(); zero_dos(); zero_hist(); zero_hist_ref();

	// setup the reference histogram
	ln_f = ln_f_initial;
	for(i = 0; i < referencesteps; i++)																					
			wanglandau();

	for(i=0; i <= modulus; i++)
			hist_ref[i] = hist[i] - reference_level;

	// iterate the Wang Landau algorithm
	newlattice(); zero_dos();
	for(ln_f = ln_f_initial; ln_f >= ln_f_final; ln_f = ln_f  / 2.0)			// proceeds over successively smaller adjustment factors
	{
		zero_hist(); n = 0;														// zero the histogram
		//printf("ln_f = %f\n",ln_f);
		do {
			if (++n > iterationlimit)
				return 0;
			for(i = 0; i < runsteps; i++)										// iterate the Wang Landau algorithm	
				wanglandau();
		} while( ! flat());														// check flatness of the histogram
	}

	normalize_dos();															// normalize the density of states

	microaverage();																// calculate the average magnetization at each energy level

	num_sums();																	// perform the numerical sums to obtain the thermodynamic quantities
	
	return -1;
}


void newlattice(void)				// Prepare a new lattice
{
	long i;
	long t = 0;

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

	stats();
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
		
void wanglandau(void)				// Select a trial cell and transition the lattice according to Wang Landau
{
	long flag = 0;
	long n, lattice_energy1, lattice_magnetiz1, t = 0;
	double p, ln_g0, ln_g1;

	n = rnd_cell();
	lattice_energy1 = lattice_energy - 2 * cell_energy(n);			// Energy of the trial configuration
	lattice_magnetiz1 = lattice_magnetiz - 2 * (lattice[n] - 1);
	ln_g0 = ln_g[indx(lattice_energy)];								// Density of states of the current and trial configurations						
	ln_g1 = ln_g[indx(lattice_energy1)];
	
	if (hist_ref[indx(lattice_energy1)] >= 0)						// Confine to reference histogram
	{
		if (ln_g1 > ln_g0)						// a flip to a configuration with a higher density of states is subject to trial
		{
				p = exp(ln_g0 - ln_g1);
				if (rnd_double() < p)
					flag = -1;
		}
		else														// a flip to a configuration with a lower density of states is always taken
			flag = -1;
	}

	if (flag)														// transition to the new configuration
	{
		flip(n);
		lattice_energy = lattice_energy1;
		lattice_magnetiz = lattice_magnetiz1;
	}
		ln_g[indx(lattice_energy)] += ln_f;
		hist[indx(lattice_energy)] += 1;
		magnetiz_avg[indx(lattice_energy)] += abs_double(lattice_magnetiz);
}

void zero_dos(void)					// initialize density of states by setting ln_g[i] = 0 (=> g[i] = 1)
{
	long i;									// project settings will the direct optimizer to inline small functions

	for (i = 0; i<modulus; i++)
	{
		ln_g[i] = 0.0;
	}
}

void microaverage(void)				// calculate the microcanonical average magnetization at each energy level
{
	long i;

	for (i = 0; i<=modulus; i++)
		if (hist_ref[i] >= 0)
			magnetiz_avg[i] /= (double) hist[i];
}


void normalize_dos(void)			// normalize the density states so that ln[ground state] = ln(2)
{
	double delta;
	long i;

	delta = ln_g[0] - ln_2;			// lowest energy state has 2 configuartions
	for (i = 0; i<=modulus; i++)
		ln_g[i] -= delta;			
}

void zero_hist(void)				// initialize histogram and the array used for the calculation of average magnetization
{
	long i;									// project settings will the direct optimizer to inline small functions

	for (i = 0; i<=modulus; i++)
	{
		hist[i] = 0;
		magnetiz_avg[i] = 0;
	}
}

void zero_hist_ref(void)			// initialize reference histogram 
{
	long i;									// project settings will the direct optimizer to inline small functions

	for (i = 0; i<=modulus; i++)
	{
		hist_ref[i] = 0;
	}
}

long flat(void)						// check the flatness of the histogram against the criterion
{
	long i, h;
	long hist_max = 0;
	long hist_min = 2147483647 ;
	double range;

	for (i = 0; i <= modulus; i++)									// scan the density of states and identify the maximum and minimum values
	{
		if (hist_ref[i] >= 0)										// ignore energy levels that cannot be reached
		{
			h = hist[i];
			if (h > hist_max)
				hist_max = h;
			if (h < hist_min)
				hist_min = h;
		}
	}

	range = ((double) (hist_max - hist_min)) / ((double) (hist_max + hist_min));				// compare the range against the flatness criterion
	//printf("range %f\tmax %d\tmin %d\n", range, hist_max, hist_min);

	if (range < (1 - flatness_criterion))
		return -1;
	else
		return 0;
}

void num_sums(void)															// perform the numerical sums over the density of states
{
	long i,j,e;
	double l, l_max;
	double e1, e2, h, h_max, f, t, sum_f;														
	double sum_e, sum_e2 = 0.0, sum_m = 0.0;
	double betastep;

	// prepare for the integration
	beta = betamin; betastep = (betamax - betamin) / (samples - 1); h_max = 0.0;
	
	// conduct an integration at each value of beta
	for (i = 0; i < samples; i++)
	{
		// pre-sum scan to remove largest factor from the exponential to improve accuracy
		l_max = 0;
		for (j = 0; j <= modulus; j++)	
			if (hist_ref[j] >= 0)												// ignore unreachable configurations
			{
				e = (4 * j - 2 * modulus);										// e = lattice energy
				l = ln_g[j] - beta * e;
				if ( l >= l_max)
					l_max = l;
			}

		// conduct the sum
		sum_f = 0.0; sum_e = 0.0; sum_e2 = 0.0; sum_m = 0.0; 
		for (j = 0; j <= modulus; j++)											
			if (hist_ref[j] >= 0)													// ignore unreachable configurations
			{
				e = (4 * j - 2 * modulus);										// e = lattice energy, E
				f = exp(ln_g[j] - beta * e - l_max);							// f = g[E]*Exp[-Beta * E]
				sum_e += e * f;
				sum_e2 += e * e * f;
				sum_m += magnetiz_avg[j] * f;
				sum_f += f;														// sum_f = partition function, Z
			}
		t = sum_f * modulus;
		e1 = sum_e / t;								// normalize with partition function and lattice size
		e2 = sum_e2 / (t * modulus); 
		h = Boltzmann * beta * beta * (e2 - e1 * e1);

		magnetiz_list[i] =  sum_m / t;
		energy_list[i] = e1;
		energy2_list[i] =e2;
		heatcapacity_list[i] = h;
		beta_list[i] = beta;

		// screening for the critical point
		if (h > h_max)
		{
			h_max = h;
			beta_critical = beta;
			sample_critical = i;
		}

		beta += betastep;
	}
}			

void stats(void)					// Calculate the stats for the lattice 
{
	long i;
	double e;
	double sum_e = 0.0;
	double sum_m = 0.0;
	double z = modulus * 2.0;		// normalization factor includes 2.0 since each interaction energy has been counted twice
	
	for(i=0; i<modulus; i++)
	{
		sum_m += lattice[i] - 1;
		e = cell_energy(i);
		sum_e += e;
	}

	// overall lattice values
	lattice_energy = (long) sum_e / 2;			// divide by 2.0 since each interaction energy has been counted twice
	lattice_magnetiz = (long) sum_m;

	// averages per site
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

void print_stats(void)			// Output the latest stats
{
	long i;

	printf("hist\tlg_g\tmag\n");
	for (i = 0; i <= modulus; i++)
		printf("%d\t%f\t%f\n",hist[i],ln_g[i],magnetiz_avg[i]);
	printf("\n");
}

void fprint_experiment(void)			// Output the results of a simulation run
{
	FILE *fp;
	int i;

	fopen_s(&fp, "wanglandau.txt","w");
	if (fp != NULL)
	{
		fprintf(fp,"Beta\t Magnetization\tHeat capacity\tEnergy\tEnergy^2\n");
		for (i = 0; i < samples; i++)
			fprintf(fp,"%f\t%f\t%f\t%f\t%f\n", beta_list[i], magnetiz_list[i], heatcapacity_list[i], energy_list[i], energy2_list[i]);
	/*for (i = 0; i < samples; i++)
		fprintf(fp,"%f\t%f\n", beta_list[i], magnetiz_list[i]);
	for (i = 0; i < samples; i++)
		fprintf(fp,"%f\t%f\n", beta_list[i], energy_list[i]);
	for (i = 0; i < samples; i++)
		fprintf(fp,"%f\t%E\n", beta_list[i], heatcapacity_list[i]);
	*/
	fclose(fp);	
	}

	fopen_s(&fp, "wanglandau3.txt","w");
	if (fp != NULL)
	{
		fprintf(fp,"Energy\tln(density of states)\n");
		for (i = 0; i <= modulus; i++)
			fprintf(fp,"%d\t%f\n", 4*i - 2*modulus, ln_g[i]);
	fclose(fp);	
	}

}

void fprint_histogram(void)				// Output the histogram of critical values of Beta
{
	FILE *fp;
	int i;

	fopen_s(&fp, "wanglandau2.txt","w");
	if (fp != NULL)
	{
		fprintf(fp,"Sample\tBeta\tCount\n");
		for (i = 0; i < samples; i++)
			fprintf(fp,"%d\t%f\t%d\n", i, beta_list[i], critical_histogram[i]);
	fclose(fp);	
	}
}