/* 
 * File:   main.cpp
 * Author: Jonathan
 *
 * Created on 03 January 2013, 00:33
 */


#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <ctime>
//#include <sys/resource.h>
using namespace std;

class measurement {
public:
    double measure;
    double temp;
};
std::vector<std::vector<int> > lattice;
std::vector<measurement> magnetisation;

// Most commonly adjusted paramaters
double TEMP_MAX = 1.0;
double TEMP_MIN = 5.0;
double NUM_STEP = 200;
int EQ_STEP = 1e7;								// runsteps * runcount
int latticesize = 32;							

double TEMP_STEP = (TEMP_MAX - TEMP_MIN)/NUM_STEP;
double TEMP = TEMP_MAX;
double CONST_J = 0.25;
int currentenergy = 0;

void initialiselattice(void);
void populatelattice(void);
void reachequilibrium(double TEMP);
void flip(int x, int y);
double sumneighbours(int x, int y);
void initialenergy(void);
void writemag(void);
void recmag(double mag);
double calcmag(void);
/*
 * 
 */
int main(int argc, char** argv) {
    srand(time(0));
    
    initialiselattice();
    cout << "Lattice Initialised..." << endl;
    populatelattice();
    cout << "Lattice Populated..." << endl;
    initialenergy();
    cout << currentenergy << endl;
    for(int i = 0; i < NUM_STEP; i++){
        reachequilibrium(TEMP);
        TEMP = TEMP - TEMP_STEP;
        cout << TEMP << endl;
    }
    writemag();
    return 0;
}

double calcmag(void){
 double mag = 0;
        for(int i = 0; i < latticesize; i++){
        for(int j = 0; j < latticesize; j++){
            mag = mag + lattice[i][j];
        }
    }
 return mag;
}
void recmag(double mag){
    measurement magnet;
    mag = mag/pow(latticesize,2.0);
    magnet.measure = mag;
    magnet.temp = TEMP;
    magnetisation.push_back(magnet);
}


void writemag(void) {
    std::ofstream output;
    output.open("mag.dat");
    for (int i = 0; i < magnetisation.size(); i++) {
        output << magnetisation[i].temp << "\t" << magnetisation[i].measure << endl;
    }
    output.close();
}

void initialenergy(void){
    currentenergy = 0;
    for(int i = 0; i < latticesize; i++){
        for(int j = 0; j < latticesize; j++){
            currentenergy = currentenergy + sumneighbours(i,j);
        }
    }
}

void flip(int x, int y){
    double energy = sumneighbours(x, y);
    if(exp((2.0*energy)/TEMP) > (double)rand()/(double)RAND_MAX){
        lattice[x][y] = -lattice[x][y];
    }
}

double sumneighbours(int x, int y) {
    double energy = 0;
    double current = lattice[x][y];
    int XMIN = 1, YMIN = 1;
    int XMAX = 2, YMAX = 2;
    /*
    if (x == 0) {
        XMIN = 0;
        if (y == 0) {
            YMIN = 0;
        } else if (y == latticesize - 1) {
            YMAX = 1;
        }
    } else if (x == latticesize - 1) {
        XMAX = 1;
        if (y == 0) {
            YMIN = 0;
        } else if (y == latticesize - 1) {
            YMAX = 1;
        }
    }

    for (int i = x - XMIN; i < x + XMAX; i++) {
        for (int j = y - YMIN; j < y + YMAX; j++) {
            energy = energy - CONST_J * current * lattice[i][j];
        }
    }
    
    
    
    return (energy + CONST_J);
*/
    if(x != 0){
        energy = energy - CONST_J*current*lattice[x - 1][y];
    }
    if(x != latticesize - 1){
        energy = energy - CONST_J*current*lattice[x + 1][y];
    }
    if(y != 0){
        energy = energy - CONST_J*current*lattice[x][y - 1];
    }
    if(y != latticesize - 1){
        energy = energy - CONST_J*current*lattice[x][y + 1];
    }
    return(energy);

}

/*
double edgecase(int x, int y){
    
}
*/

void reachequilibrium(double TEMP){
    double mag = 0;
    int x, y;
    //double energy;
    int j = 0;
    double nummeasures = 0;
    for (int i = 0; i < EQ_STEP; i++){
        x = (int)(rand() % latticesize);
        y = (int)(rand() % latticesize);
        flip(x,y);
        j++;
        if( j == 1e3){						// HERE		runsteps
            j = 0;
            mag += abs(calcmag());
            nummeasures++;
        }
    }
    mag = mag/nummeasures;
    recmag(mag);
}

void initialiselattice(void){
    lattice.resize(latticesize);
    for(int i = 0; i < latticesize; i ++){
        lattice[i].resize(latticesize);
    }
}

void populatelattice(void){
    for(int i = 0; i < latticesize; i++){
        for(int j = 0; j < latticesize; j++){
            if((double)rand()/(double)RAND_MAX > 0.5){
                lattice[i][j] = 1;
            }
            else{
                lattice[i][j] = -1;
            }
        }
    }
}