#include <math.h>
#include "vectors.h"

// Function to calculate the norm of a vector
double vectorNorm(const Vector3D* v) {
    return sqrt(v->x * v->x + v->y * v->y + v->z * v->z);
}

// Function to calculate the dot product of two vectors
double dotProduct(const Vector3D* v1, const Vector3D* v2) {
    return v1->x * v2->x + v1->y * v2->y + v1->z * v2->z;
}

void addVectors(const Vector3D* v1, const Vector3D* v2, Vector3D* result) {
  result->x = v1->x + v2->x;
  result->y = v1->y + v2->y;
  result->z = v1->z + v2->z;
}

void copyVectors(const Vector3D* v, Vector3D* result) {
  result->x = v->x;
  result->y = v->y;
  result->z = v->z;
}

// Function to calculate the cross product of two vectors
void crossProduct(const Vector3D* v1, const Vector3D* v2, Vector3D* result) {
    result->x = v1->y * v2->z - v1->z * v2->y;
    result->y = v1->z * v2->x - v1->x * v2->z;
    result->z = v1->x * v2->y - v1->y * v2->x;
}

// Function to perform scalar-vector multiplication
void scalarVectorMultiply(double scalar, const Vector3D* v, Vector3D* result) {
    result->x = scalar * v->x;
    result->y = scalar * v->y;
    result->z = scalar * v->z;
}

// Function to rotate a vector around a given axis by a specified angle (in radians)
void rotateVector(const Vector3D* v, const Vector3D* axis, double angle, Vector3D* result) {
    double cosTheta = cos(angle);
    double sinTheta = sin(angle);
    double dot = dotProduct(axis, v); // Compute dot product of axis and vector

    // Rodrigues' rotation formula
    result->x = cosTheta * v->x + sinTheta * (axis->y * v->z - axis->z * v->y) + (1 - cosTheta) * dot * axis->x;
    result->y = cosTheta * v->y + sinTheta * (axis->z * v->x - axis->x * v->z) + (1 - cosTheta) * dot * axis->y;
    result->z = cosTheta * v->z + sinTheta * (axis->x * v->y - axis->y * v->x) + (1 - cosTheta) * dot * axis->z;
}

// Function to calculate the distance between two vectors
double distanceBetweenVectors(const Vector3D* v1, const Vector3D* v2) {
    double dx = v1->x - v2->x;
    double dy = v1->y - v2->y;
    double dz = v1->z - v2->z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}


// Function to calculate the arc length on a sphere of radius R between two points v1 and v2
double arcLengthOnSphere(double R, const Vector3D* v1, const Vector3D* v2) {
    // Convert Cartesian coordinates to spherical coordinates (latitude and longitude)
    double lat1 = atan2(v1->z, sqrt(v1->x * v1->x + v1->y * v1->y));
    double lon1 = atan2(v1->y, v1->x);
    double lat2 = atan2(v2->z, sqrt(v2->x * v2->x + v2->y * v2->y));
    double lon2 = atan2(v2->y, v2->x);

    double deltaLon = lon2 - lon1;
    double centralAngle = acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(deltaLon));
    return R * centralAngle;
}
