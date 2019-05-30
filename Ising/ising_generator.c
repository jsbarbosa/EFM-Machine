#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "headers.h"

int **spins;

int main(int argc, char const *argv[])
{
  // sampleTemperatures("results.dat", 1e-6, 1e1, 5);

  // double E, mag;
  // createSpins();
  // sampleTemperature(1e-6, &E, &mag);
  // clean();

  doEvolution(1e-6, 1e6, 1e3);
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
  double E, mag;
  double t;

  for(t = t_min; t < t_max; t += (t_max - t_min) / n)
  {
    createSpins();
    sampleTemperature(t, &E, &mag);
    fprintf(file, "%e %e %e\n", t, E, mag);
    clean();
  }
  fclose(file);
}

void sampleTemperature(double T, double *E, double *mag)
{
  int i;
  *E = 0;
  *mag = 0;
  for(i = 0; i < N_ignore; i++) evolve(T);
  for(i = 0; i < N_samples; i++)
  {
    evolve(T);
    *E += getEnergy();
    *mag += getMagnitude();
  }
  *E = *E / N_samples;
  *mag = *mag / N_samples;
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
  double center, up, down, left, right;
  center = energyAt(i, j);
  up = energyAt(posBoundaries(i + 1, Ly), j);
  down = energyAt(posBoundaries(i - 1, Ly), j);
  left = energyAt(i, posBoundaries(j - 1, Lx));
  right = energyAt(i, posBoundaries(j + 1, Lx));

  return center + up + down + left + right;
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
