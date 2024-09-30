#ifndef VECTORS_H
#define VECTORS_H

typedef struct {
    double x;
    double y;
    double z;
} Vector3D;

typedef struct {
    double theta;
    double phi;
} Vector2DSphere;

double vectorNorm(const Vector3D* v);
double dotProduct(const Vector3D* v1, const Vector3D* v2);
void addVectors(const Vector3D* v1, const Vector3D* v2, Vector3D* res);
void copyVectors(const Vector3D* v, Vector3D* res);
void crossProduct(const Vector3D* v1, const Vector3D* v2, Vector3D* res);
void scalarVectorMultiply(double scalar, const Vector3D* v, Vector3D* res);
void rotateVector(const Vector3D* v, const Vector3D* axis, double angle, Vector3D* res);
double distanceBetweenVectors(const Vector3D* v1, const Vector3D* v2);
double distanceBetweenVectorsPolar(const double r, const Vector2DSphere* v1, const Vector2DSphere* v2);
double arcLengthOnSphere(double R, const Vector3D* v1, const Vector3D* v2);

#endif // VECTORS_H

