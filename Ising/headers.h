void createSpins(void);
void printSpins(void);
double energyAt(int, int);
int posBoundaries(int, int);
void evolve(double);
double associatedEnergy(int, int);
double normalRandom(void);
double getEnergy(void);
void sampleTemperature(double T, double *, double *, double *, double *);
int getMagnitude(void);

void sampleTemperatures(const char *filename, double t_min, double t_max, int n);
void clean(void);
void printFile(const char *filename);
void doEvolution(double T, int n, int);

#define Lx 8
#define Ly 8
#define J 1

#define N_ignore 1e2
#define N_samples 1e4

#define T_MIN 1e-6
#define T_MAX 10
#define N_temperatures 5
