void createSpins(void);
void printSpins(void);
double energyAt(int, int);
int posBoundaries(int, int);
void evolve(double);
double associatedEnergy(int, int);
double normalRandom(void);
double getEnergy(void);
void sampleTemperature(double T, double *E, double *mag);
int getMagnitude(void);

void sampleTemperatures(const char *filename, double t_min, double t_max, int n);
void clean(void);
void printFile(const char *filename);
void doEvolution(double T, int n, int);

#define Lx 200
#define Ly 200
#define J 1

#define N_ignore 1e3
#define N_samples 1e5

#define T_MIN 1e-6
#define T_MAX 10
#define N_temperatures 5
