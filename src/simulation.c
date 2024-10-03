#include "simulation.h"
#include "config.h"
#include "vectors.h"

#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void diff(double* arr, double* result, int length) {
  for (int i = 0; i < length - 1; i++) {
      result[i] = arr[i + 1] - arr[i];
  }
}

void compForceList(double* rr,
                   double* Vpot, double* rF,
                   double* Flist, int lenght) {

  double* diffRR = (double*)malloc((lenght-1) * sizeof(double));
  double* diffVpot = (double*)malloc((lenght-1) * sizeof(double));

  // Calculate differences
  diff(rr, diffRR, lenght);
  diff(Vpot, diffVpot, lenght);

  // Calculate r_F
  for (int i = 0; i < lenght - 1; i++) {
      rF[i] = diffRR[i] / 2.0 + rr[i];
  }

  // Calculate F_list
  for (int i = 0; i < lenght - 1; i++) {
      Flist[i] = -diffVpot[i] * 1.60218e-19 / (diffRR[i] * 1e-10)*6.022e26*1e10*1e-30;
  }

  // Free allocated memory
  free(diffRR);
  free(diffVpot);
}

void initRandomParticle(Particle *particles,
                        int       noP,
                        double radius,
                        double vel) {
  // particles is a array of 2 * noP particles.

  srand (time(NULL));
  const double pi = 4 * atan(1);
  double distance = 0.0;
  Particle *particle1 = NULL;
  Particle *particle2 = NULL;
  double theta = 0.0;
  double phi = 0.0;

  for (int i = 0; i < noP; ++i) {
    particle1 = particles + 2*i;
    particle2 = particles + 2*i + 1;

    particle1->index = i;
    particle2->index = i;
    particle1->time = 0;
    particle2->time = 0;
    particle1->success = 0;
    particle2->success = 0;

    // Generate two random positions on a sphere
    // such that distance between two positions is
    // greater than ICD Distance.

    do {
      theta = pi * ((double)rand() / RAND_MAX);
      phi = 2 * pi * ((double)rand() / RAND_MAX);

      particle1->pos.x = radius * sin(theta) * cos(phi);
      particle1->pos.y = radius * sin(theta) * sin(phi);
      particle1->pos.z = radius * cos(theta);

      theta = pi * ((double)rand() / RAND_MAX);
      phi = 2 * pi * ((double)rand() / RAND_MAX);
      particle2->pos.x = radius * sin(theta) * cos(phi);
      particle2->pos.y = radius * sin(theta) * sin(phi);
      particle2->pos.z = radius * cos(theta);

      distance = distanceBetweenVectors (&(particle1->pos), &(particle2->pos));

    } while( distance < icddist);

    // Generate angular velocity of particle 1
    Vector3D ranVector = {0, 1, 0};
    Vector3D postAng;

    crossProduct(&ranVector, &(particle1->pos), &postAng);

    scalarVectorMultiply((1.0 / vectorNorm(&postAng)) * (vel / radius),
                         &postAng, &(particle1->angVel));

    // Generate angular velocity of particle 2

    crossProduct(&ranVector, &(particle2->pos), &postAng);

    scalarVectorMultiply((1.0 / vectorNorm(&postAng)) * (vel / radius),
                         &postAng, &(particle2->angVel));
                         
    particle1->torque.x = 0;
    particle1->torque.y = 0;
    particle1->torque.z = 0;

    particle2->torque.x = 0;
    particle2->torque.y = 0;
    particle2->torque.z = 0;
  }

}


void moveParticle(Particle *particle, Vector3D *angVel) {
    double rotAngle = 0.0;
    Vector3D meanAngVel = {0., 0., 0.};
    Vector3D rotVector = {0., 0., 0.};

    // Move particle according to angular velocity
    addVectors (&(particle->angVel), angVel, &meanAngVel);
    scalarVectorMultiply (0.5, &meanAngVel, &meanAngVel);
    rotAngle = vectorNorm (&meanAngVel) * dt;
    scalarVectorMultiply (dt/rotAngle, &meanAngVel, &rotVector);

    // Rotate to move particle
    rotateVector (&(particle->pos), &rotVector, rotAngle, &(particle->pos));
    copyVectors (&(particle->angVel), angVel);
}

void updateParticle(Particle *particle, Vector3D* force,
                    Vector3D *oldAngVel,
                    Vector3D *oldDelAngVel,
                    double radius) {

    double epsilon = 1e-9;
    Vector3D delAngVel = {0, 0, 0};
    Vector3D meanDelAngVel = {0., 0., 0.};

    // Calculate Angular velocity of particle
    // under force
    if (vectorNorm (force) > epsilon) {
      // Calculate torque on particle
      crossProduct (&(particle->pos), force, &(particle->torque));

      // Change in angular velocity of particle
      scalarVectorMultiply(dt / (mHe*radius*radius),
                            &(particle->torque), &delAngVel);

      if (vectorNorm (oldDelAngVel) > epsilon) {
        addVectors (oldDelAngVel, &delAngVel , &meanDelAngVel);
        scalarVectorMultiply (0.5*dt, &meanDelAngVel, &meanDelAngVel);
      }
      else {
        scalarVectorMultiply (dt, &delAngVel, &meanDelAngVel);
      }

      addVectors (oldAngVel, &meanDelAngVel, &(particle->angVel));

      // Reset old values
      copyVectors (&delAngVel, oldDelAngVel);
    }

    if (vectorNorm (&(particle->angVel)) > epsilon) moveParticle(particle, oldAngVel);
}


int testWithoutForce(Particle *particle1,
                     Particle *particle2)  {
  int controlpar = 0;
  int time = 0;
  int dtTest = 1000;

  double controldistance = 20.0;
  double distance = 0.0;

  Particle* copy1 = (Particle *)malloc (sizeof (Particle));
  Particle* copy2 = (Particle *)malloc (sizeof (Particle));

  // Copy old particles
  memcpy(copy1, particle1, sizeof(Particle));
  memcpy(copy2, particle2, sizeof(Particle));

  while (time < maxtime && !controlpar) {
    distance = distanceBetweenVectors(&(copy1->pos), &(copy2->pos));
    if (distance < controldistance) controlpar = 1;

    // Update particle 1
    moveParticle(copy1, &(copy1->angVel));

    // Update particle 2
    moveParticle(copy2, &(copy2->angVel));

    time += dtTest;
  }

  return controlpar;

}

void particleSimulation(int index,
                        Particle *particle1,
                        Particle *particle2,
                        double radius,
                        double *Flist) {

  int rkFlag = 0;
  int saveflag = 0;
  int time = 0;

  Vector3D oldAngVel1= {0., 0., 0.};
  Vector3D oldAngVel2= {0., 0., 0.};
  Vector3D oldDelAngVel1 = {0, 0, 0};
  Vector3D oldDelAngVel2 = {0, 0, 0};

  double indxd = 0.0;
  int indx = 0;
  double distance = 0.0;
  double force = 0.0;
  Vector3D forc1 = {0, 0, 0};
  Vector3D forc2 = {0, 0, 0};
  
  char affname1[80];
  char affname2[80];
  FILE *ffile1 = NULL;
  FILE *ffile2 = NULL;
  
  if (saveflag) {
    sprintf (affname1, "./results/pair_%d_particle1.txt", index);
    ffile1 = fopen(affname1, "w");

    sprintf (affname2, "./results/pair_%d_particle2.txt", index);
    ffile2 = fopen(affname2, "w");
  }

  // Check if particle meets without using force
  if(!testWithoutForce(particle1, particle2))  {
    while(time < maxtime) {

      if (saveflag) {
        // Write particle properties after the simulation
        fprintf (ffile1,"%d\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n",
                 particle1->success, particle1->time,
                 particle1->pos.x, particle1->pos.y, particle1->pos.z);
        fprintf (ffile2,"%d\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n",
                 particle1->success, particle1->time,
                 particle2->pos.x, particle2->pos.y, particle2->pos.z);
      }


      if (!rkFlag) {
        // Determine the force between both particles
        distance = distanceBetweenVectors (&(particle1->pos), &(particle2->pos));
        indxd = (distance - 2.005) / 0.01;
        if (indxd < 22.0) indx = (int)(indxd);
        else indx = 2000;

        if (indx > 1799) indx = 1799;
        force = Flist[indx];

        // Force on particle 1
        forc1.x = force * (particle1->pos.x - particle2->pos.x) / distance;
        forc1.y = force * (particle1->pos.y - particle2->pos.y) / distance;
        forc1.z = force * (particle1->pos.z - particle2->pos.z) / distance;

        // Force on particle 2
        scalarVectorMultiply (-1, &forc1, &forc2);

        // Update Particle 1
        updateParticle (particle1, &forc1, &oldAngVel1, &oldDelAngVel1, radius);

        // Update Particle 2
        updateParticle (particle2, &forc2, &oldAngVel2, &oldDelAngVel2, radius);

        if ( distanceBetweenVectors (&(particle1->pos), &(particle2->pos)) < icddist) {
          particle1->success = 1;
          particle2->success = 1;
          break;
        }
      }
      time += dt;
      particle1->time = time;
      particle2->time = time;

    }
  }
  
  if (saveflag) {
    fclose(ffile1);
    fclose(ffile2);
  }

}
