#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include "simulation.h"
#include "config.h"

void printUsage(const char *progName) {
  fprintf(stderr, "Usage: %s -h <he_number> -v <velocity> -n <number>\n", progName);
  exit(EXIT_FAILURE);
}

int main(int argc, char *argv[]) {
  int opt;
  int helium_number = -1;
  double radius = -1.0;
  double velocity = -1.0;
  int number = -1;

  while ((opt = getopt(argc, argv, "h:v:n:")) != -1) {
    switch (opt) {
      case 'h':
        helium_number = atof(optarg);
        break;
      case 'v':
        velocity = atof(optarg);
        break;
      case 'n':
        number = atoi(optarg);
        break;
      default: /* '?' */
        printUsage(argv[0]);
    }
  }

  // Check if both arguments were provided
  if (helium_number == -1 || velocity == -1.0 || number == -1) {
    printUsage(argv[0]);
  }


  char fitfilename[80] = "./data/fit_potential.txt";
  FILE *fitfile = fopen(fitfilename, "r");

  if (fitfile == NULL) {
      perror("Error opening file");
      return 1;
  }

  // Read the potential file and store the values into arrays
  double rr[1800]; // radius array
  double yss[1800]; // Singlet_Single array
  double yst[1800]; // Singlet_Triplet array
  double ytt[1800]; // Triplet_Triplet array

  int i = 0;
  while (fscanf(fitfile, "%lf\t%lf\t%lf\t%lf\n", &rr[i], &yss[i], &yst[i], &ytt[i]) == 4) {
      i++;
  }

  fclose(fitfile);

  // Compute interatomic force list
  double rF[1799];
  double Flist[1799];
  compForceList (rr, yss, rF, Flist, 1800);
  
  // Save forcelist
  char fffname[80] = "./results/forcelist.txt";
  FILE *fffile = fopen(fffname, "w");
  fprintf(fffile, "Radius\tForce\n");
  for (int i = 0; i < 1800; ++i) {
    fprintf(fffile, "%lf\t%lf\n", rF[i], Flist[i]);
  }

  // Array to hold twin particles
  int noP = number; // No of particle pairs
  radius = 2.22 * pow((double)helium_number, 1.0/3.0);
  Particle* particles = (Particle *)malloc (2 * noP * sizeof (Particle));

  // Initialize random particleSimulation
  initRandomParticle (particles, noP, radius, velocity);

  // Write particle properties in a file
  char ffname[80];
  sprintf(ffname, "./results/particlefile.txt");
  FILE *ffile = fopen(ffname, "w");
  fprintf (ffile,"Index\tSuccess\tTime\tX\tY\tZ\twx\twy\twz\ttx\tty\ttz\n");
  for (int i = 0; i < 2*noP; ++i) {
    fprintf (ffile,"%d\t%d\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n",
             particles[i].index, particles[i].success, particles[i].time,
             particles[i].pos.x, particles[i].pos.y, particles[i].pos.z,
             particles[i].angVel.x, particles[i].angVel.y, particles[i].angVel.z,
             particles[i].torque.x, particles[i].torque.y, particles[i].torque.z);
  }
  fclose(ffile);

  // Simulate particles
  #pragma omp parallel for
  for (int i = 0; i < noP; ++i) {
    particleSimulation (i, particles + 2*i, particles + 2*i + 1,
                        radius, Flist);
  }

  // Write particle properties after the simulation
  char affname[80];
  sprintf(affname, "./results/particlefileafter.txt");
  ffile = fopen(affname, "w");
  fprintf (ffile,"Index\tSuccess\tTime\tX\tY\tZ\twx\twy\twz\ttx\tty\ttz\n");
  for (int i = 0; i < 2*noP; ++i) {
    fprintf (ffile,"%d\t%d\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n",
             particles[i].index, particles[i].success, particles[i].time,
             particles[i].pos.x, particles[i].pos.y, particles[i].pos.z,
             particles[i].angVel.x, particles[i].angVel.y, particles[i].angVel.z,
             particles[i].torque.x, particles[i].torque.y, particles[i].torque.z);
  }
  fclose(ffile);

  free(particles);
  return 0;
}
