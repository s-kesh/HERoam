#ifndef SIMULATION_H
#define SIMULATION_H

#include "vectors.h"

typedef struct {
  int index;
  int success;
  double time;
  Vector3D pos;
  Vector3D angVel;
  Vector3D torque;
} Particle;

double calculatePotential(double *pot_array,
                          double  radius);
void diff(double* arr, double* result, int length);
void compForceList(double* rr,
                   double* Vpot, double* rF,
                   double* Flist, int lenght);
void initRandomParticle(Particle *particle,
                        int       noP,
                        double radius,
                        double vel);

void moveParticle(Particle *particle, Vector3D *oldAngVel);

void updateParticle(Particle *particle, Vector3D* force,
                    Vector3D *oldAngVel,
                    Vector3D *oldDelAngVel,
                    double radius);

int testWithoutForce(Particle *particle1,
                      Particle *particle2);

void particleSimulation(int index, Particle *particle1, Particle *particle2,
                        double radius, double *Flist);

#endif
