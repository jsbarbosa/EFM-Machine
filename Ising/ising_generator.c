#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "headers.h"

int **spins;

const double TC = 2.0 / log(1.0 + sqrt(2)) * J;

int main(int argc, char const *argv[])
{
  sampleTemperatures("results.dat", 0.1, 7, 150);

  // double E, mag;
  // createSpins();
  // sampleTemperature(1e-6, &E, &mag);
  // clean();

  // doEvolution(1e-6, 1e6, 1e3);
  return 0;
}

void doEvolution(double T, int n, int print_every)
{
  int i, j = 0;
  createSpins();
  char buffer[50];

  for(i = 0; i < n; i++)
  {
    evolve(T);
    if(i % print_every == 0)
    {
      sprintf(buffer, "results/%d.ising", j);
      printFile(buffer);
      j += 1;
    }
  }
  clean();
}

void printFile(const char *filename)
{
  int i, j;
  FILE *file = fopen(filename, "w");
  for(i = 0; i < Ly; i++)
  {
    for(j = 0; j < Lx; j++)
    {
      fprintf(file, "%d ", spins[i][j]);
    }
    fprintf(file, "\n");
  }
  fclose(file);
}

void sampleTemperatures(const char *filename, double t_min, double t_max, int n)
{
  FILE *file = fopen(filename, "w");
  double E, mag, cv, chi;
  double t;
  createSpins();
  for(t = t_min; t < t_max; t += (t_max - t_min) / n)
  {
    // createSpins();
    sampleTemperature(t, &E, &mag, &cv, &chi);
    // clean();
    fprintf(file, "%e %e %e %e %e \n", t / TC, E, mag, cv, chi);
  }
  fclose(file);
  clean();
}

void sampleTemperature(double T, double *E, double *mag, double *cv, double *chi)
{
  int i;
  *E = 0;
  *mag = 0;

  int m = Lx * Ly;
  double E2 = 0, mag2 = 0;

  for(i = 0; i < N_ignore; i++) evolve(T);
  for(i = 0; i < N_samples; i++)
  {
    evolve(T);
    *E += getEnergy() / (N_samples * m);
    E2 += pow(*E, 2);
    *mag += fabs(getMagnitude()) / (N_samples * m);;
    mag2 += pow(*mag, 2);
  }
  *cv = (E2 - pow(*E, 2)) / (T * T);
  *chi = (mag2 - pow(*mag, 2)) / T;
}

void evolve(double T)
{
  int i = rand() % Ly;
  int j = rand() % Lx;
  double E = associatedEnergy(i, j);
  spins[i][j] *= -1; // change spin
  E = associatedEnergy(i, j) - E; //final - initial
  if((E >= 0) && (rand() / (double) RAND_MAX >= exp(-E / T)))
  {
    spins[i][j] *= -1;
  }
}

double energyAt(int i, int j)
{
  int k;
  int sum = 0;
  for(k = -1; k < 2; k += 2)
  {
    sum += spins[posBoundaries(i + k, Ly)][j] + spins[i][posBoundaries(j + k, Lx)];
  }
  return -J * sum * spins[i][j];
}

int getMagnitude(void)
{
  int i, j, mag = 0;
  for(i = 0; i < Ly; i++)
  {
    for(j = 0; j < Lx; j++)
    {
      mag += spins[i][j];
    }
  }
  return mag;
}

double getEnergy(void)
{
  int i, j;
  double E = 0;
  for(i = 0; i < Ly; i++)
  {
    for(j = 0; j < Lx; j++)
    {
      E += energyAt(i, j);
    }
  }
  return E;
}

int posBoundaries(int i, int size)
{
  int pos = i % size;
  if(pos < 0) pos += size;
  return pos;
}

double associatedEnergy(int i, int j)
{
  double center = energyAt(i, j);
  // double up, down, left, right;
  // up = energyAt(posBoundaries(i + 1, Ly), j);
  // down = energyAt(posBoundaries(i - 1, Ly), j);
  // left = energyAt(i, posBoundaries(j - 1, Lx));
  // right = energyAt(i, posBoundaries(j + 1, Lx));

  return center;
  // return center + up + down + left + right;
}

void createSpins(void)
{
  int i, j;
  time_t t;
  srand((unsigned) time(&t));
  spins = malloc(Ly * sizeof(int *));
  for(i = 0; i < Ly; i++)
  {
    spins[i] = malloc(Lx * sizeof(int));
  }

  for(i = 0; i < Ly; i++)
  {
    for(j = 0; j < Lx; j++)
    {
      spins[i][j] = 2 * (rand() % 2) - 1;
    }
  }
}

void printSpins(void)
{
  int i, j;
  for(i = 0; i < Ly; i++)
  {
    for(j = 0; j < Lx; j++)
    {
      printf("%d ", spins[i][j]);
    }
    printf("\n");
  }
}

void clean(void)
{
  int i;
  for(i = 0; i < Ly; i++)
  {
    free(spins[i]);
  }
  free(spins);
}
