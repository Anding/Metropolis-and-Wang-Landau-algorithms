#include <stdio.h>
#include <stdlib.h>

// Helpful resource
// http://msdn.microsoft.com/en-us/library/fw5abdx6.aspx

int dimension;
int modulus;
double Boltzmann;
char lattice[16384];

unsigned long rnd(void);			// Return a random integer between 0 and 2^32-1
void newlattice(int n);	// Prepare an n*n lattice - n MUST BE A POWER OF 2
void render(void);		// Display the lattice

int main()
{
	int i;
	for (i = 0; i<1; i++)
	{
		newlattice(16);
		render();
	}
	getchar();
	return 0;
}

unsigned long rnd()
{
	static unsigned long X;
	return X = (unsigned long) ((unsigned long long)X * 3141592621ULL + 1ULL);
	// use the 64 integer type ULL to avoid premature truncation
	// the modulus is provided by the typecast (unsigned long)
}

void newlattice (int n)
{
	int i;

	dimension = n;
	modulus = n*n;

	for (i=0; i<=modulus; ++i)
	{
		if (rnd() & (1<<31))
			lattice[i] = 2;
		else
			lattice[i] = 0;
	}

	Boltzmann = 0;
}

void render(void)
{
	int i, j, n;
	n = 0;

	for (j=0; j<dimension; ++j)
	{
		for (i=0; i<dimension; ++i)
		{
			if (lattice[n++] == 2)
				printf("+ ");
			else
				printf("- ");
		}
		printf("\n");
	}
}