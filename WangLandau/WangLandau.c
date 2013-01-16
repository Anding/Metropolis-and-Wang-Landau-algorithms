#include <stdio.h>
#include <math.h>

// Most commonly adjusted paramaters
#define dimensions 2				// maximum 4, assuming a span of 64
#define span 16						// width of lattice along every dimension
#define Tmin 1.0					// minimum value of beta
#define Tmax 10.0					// maximum value of beta
#define samples 200					// # of samples between betamin and betamax
#define	experiments 1				// # of separate experiments to compile
#define ln_f_limit 0.00001			// limit value of the adjustment factor
#define	referencesteps 10000000		// number of Monte Carlo steps to establish the reference histogram

// Other paramaters and constants
#define coldstart 0					// -1 = coldstart, 0 = hotstart
#define reference_level 100			// minimum count in the reference histogram to include an energy level
#define ln_f_initial 1.0			// initial value of the adjustment factor
#define runsteps 100000				// number of Monte Carlo steps between each check of the histogram
#define flatness_criterion 0.80		// criterion for testing the flatness of the histogram
#define iterationlimit 1000000		// criterion for avoiding stuck random walks
#define Boltzmann 1					// Boltzmann constant (used only for scaling heat capacity)
#define ln_2 0.6931471806			// ln_g[lowest energy configuration] = ln(2), for normalization

// Macro for addressing the density of states array with an energy level
#define indx(E)	(E + dimensions * modulus) / 4
#define unindx(i) (4 * i - dimensions * modulus)
	
// Global variables
double beta;						// Thermodynamic beta  = 1 / T
//double beta_critical;				// Beta at the critical point
long sample_critical;				// sample number of critical beta
long modulus;						// number of cells in the lattice
long energy_levels;					// number of energy levels
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
double Temp_list[samples];
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

int main()
{
	long i;
	unsigned long x = rnd();

	for (i = 0; i < samples; i++)				// zero the histogram
	{
		critical_histogram[i] = 0;
	}

	for (i = 0; i < experiments; i++)			// perform the experiment set
	{
		printf("%d...",i);
		while (x = rnd(), ! experiment())
			printf("\nHang with seed %u\n", x);

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
	energy_levels = (dimensions * modulus) / 2 + 1;

	// setup the reference histogram
	ln_f = ln_f_initial;
	while (hist[0] <= reference_level)											// ensure that the ground state gets on the list!
		for(i = 0; i < referencesteps; i++)																					
				wanglandau();

	for(i=0; i < energy_levels; i++)
			hist_ref[i] = hist[i] - reference_level;

	// initialize the lattice to an allowed configuration
	zero_dos();
	do {
		newlattice();
		stats();
	}
	while (hist_ref[indx(lattice_energy)] < 0);

	// iterate the Wang Landau algorithm
	for(ln_f = ln_f_initial; ln_f >= ln_f_limit; ln_f = ln_f  / 2.0)			// proceeds over successively smaller adjustment factors
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

long flat(void)						// check the flatness of the histogram against the criterion
{
	long i, h;
	long hist_max = 0;
	long hist_min = 2147483647 ;
	double range;

	for (i = 0; i < energy_levels; i++)									// scan the density of states and identify the maximum and minimum values
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
	double Temp, Tstep;

	// prepare for the integration
	Tstep = (Tmax - Tmin) / (samples - 1); Temp = Tmin;  h_max = 0.0; 
	
	// conduct an integration at each value of beta
	for (i = 0; i < samples; i++)
	{
		// set beta
		beta = 1.0 / Temp;

		// pre-sum scan to remove largest factor from the exponential to improve accuracy
		l_max = 0;
		for (j = 0; j < energy_levels; j++)	
			if (hist_ref[j] >= 0)												// ignore unreachable configurations
			{
				e = unindx(j);													// e = lattice energy									
				l = ln_g[j] - beta * e;
				if ( l >= l_max)
					l_max = l;
			}

		// conduct the sum
		sum_f = 0.0; sum_e = 0.0; sum_e2 = 0.0; sum_m = 0.0; 
		for (j = 0; j < energy_levels; j++)											
			if (hist_ref[j] >= 0)													// ignore unreachable configurations
			{
				e = unindx(j);													// e = lattice energy, E									
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
		Temp_list[i] = Temp;

		// screening for the critical point
		if (h > h_max)
		{
			h_max = h;
			sample_critical = i;
		}

		Temp += Tstep;
	}
}		

void zero_dos(void)					// initialize density of states by setting ln_g[i] = 0 (=> g[i] = 1)
{
	long i;									// project settings will the direct optimizer to inline small functions

	for (i = 0; i < energy_levels; i++)
	{
		ln_g[i] = 0.0;
	}
}

void zero_hist(void)				// initialize histogram and the array used for the calculation of average magnetization
{
	long i;									// project settings will the direct optimizer to inline small functions

	for (i = 0; i < energy_levels; i++)
	{
		hist[i] = 0;
		magnetiz_avg[i] = 0;
	}
}

void zero_hist_ref(void)			// initialize reference histogram 
{
	long i;									// project settings will the direct optimizer to inline small functions

	for (i = 0; i < energy_levels; i++)
	{
		hist_ref[i] = 0;
	}
}

void normalize_dos(void)			// normalize the density states so that ln[ground state] = ln(2)
{
	double delta;
	long i;

	delta = ln_g[0] - ln_2;			// lowest energy state has 2 configuartions
	for (i = 0; i < energy_levels; i++)
		ln_g[i] -= delta;			
}

void microaverage(void)				// calculate the microcanonical average magnetization at each energy level
{
	long i;

	for (i = 0; i < energy_levels; i++)
		if (hist_ref[i] >= 0)
			magnetiz_avg[i] /= (double) hist[i];
}

void print_stats(void)			// Output the latest stats
{
	long i;

	printf("E\thist\tlg_g\tmag\n");
	for (i = 0; i < energy_levels; i++)
		printf("%d\t%d\t%f\t%f\n",unindx(i), hist[i], ln_g[i], magnetiz_avg[i]);
	printf("\n");
}

void fprint_experiment(void)			// Output the results of a simulation run
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

	fopen_s(&fp, "density.txt","w");
	if (fp != NULL)
	{
		fprintf(fp,"Energy\tln(density of states)\n");
		for (i = 0; i < energy_levels; i++)
			fprintf(fp,"%d\t%f\t%f\n", unindx(i), ln_g[i], magnetiz_avg[i]);
		fclose(fp);	
	}

	fopen_s(&fp, "signature.txt","w");
	if (fp != NULL)
	{
		fprintf(fp,"Dimensions = %d\n",dimensions);
		fprintf(fp,"Span = %d\n",span);
		fprintf(fp,"Reference histogram steps = %d\n",referencesteps);
		fprintf(fp,"ln f_limit = %f\n",ln_f_limit);
		fprintf(fp,"No. experiments = %d\n",experiments);
		fclose(fp);	
	}

}

void fprint_histogram(void)				// Output the histogram of critical values of Beta
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

// Library routines

// Return a random integer x between 0 and 1<<31
unsigned long rnd(void)				
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

// Return a random integer 0 <= x < n
unsigned long rnd_long(long n)		
{
	return (unsigned long) (rnd_double() * (double) n);
}

// abs for double precision numbers
double abs_double(double x)			
{
	if (x > 0)
		return x;
	else
		return -x;
}

// Ising routines

// Prepare a new lattice
void newlattice(void)				
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

// Pick a random cell in the lattice
long rnd_cell(void)					
{
	return rnd_long(modulus);
}

// Find the nearest neighbours of the given cell
void find_neighbours(long n)		
{
	int i;

	for (i = 0; i<dimensions; i++)
	{
		// in each dimension, i, there are two neighbours at a distance reach[i] away
		neighbours[2*i] = (n + reach[i] + modulus) % modulus;		// make sure to bring the offset back into range
		neighbours[2*i+1] = (n - reach[i] + modulus) % modulus;
	}
}

// Sum the states of the neighbours
long sum_neighbours(void)			
{
	long i;
	long e=0;

	for (i=0; i < 2*dimensions; i++)
		e += lattice[neighbours[i]] - 1;

	return e;
}

// Sum the interaction energy attributable to cell n
long cell_energy(long n)			
{
	find_neighbours(n);
	return -1 *(lattice[n] - 1) * sum_neighbours();
}

// Flip cell n
void flip(long n)					
{
	lattice[n] = 2 - lattice[n];
}

// Calculate the stats for the lattice 
void stats(void)					
{
	long i;
	double e;
	double sum_e = 0.0;
	double sum_m = 0.0;
	double z = modulus * 2.0;					// normalization factor includes 2.0 since each interaction energy has been counted twice
	
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

// Display the lattice
void render(void)					
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

// Helpful C resource
// http://msdn.microsoft.com/en-us/library/fw5abdx6.aspx
