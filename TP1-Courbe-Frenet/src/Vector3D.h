#ifndef VECTOR3D_H
#define VECTOR3D_H

#include <cmath>

// Structure pour représenter un point 3D
struct Vector3D {
    double x, y, z;

    // Opérateur d'addition
    Vector3D operator+(const Vector3D& other) const {
        return Vector3D{ x + other.x, y + other.y, z + other.z };
    }

    // Opérateur de soustraction
    Vector3D operator-(const Vector3D& other) const {
        return Vector3D{ x - other.x, y - other.y, z - other.z };
    }

    // Opérateur de multiplication par un scalaire
    Vector3D operator*(double scalar) const {
        return Vector3D{ x * scalar, y * scalar, z * scalar };
    }

    // Opérateur de division par un scalaire
    Vector3D operator/(double scalar) const {
        return Vector3D{ x / scalar, y / scalar, z / scalar };
    }

    // Calcul de la norme
    double norm() const {
        return sqrt(x * x + y * y + z * z);
    }

    // Normalisation du vecteur
    Vector3D normalize() const {
        double n = norm();
        if (n == 0) return Vector3D{0, 0, 0};
        return (*this) / n;
    }

    // Produit vectoriel
    Vector3D cross(const Vector3D& other) const {
        return Vector3D{
            y * other.z - z * other.y,
            z * other.x - x * other.z,
            x * other.y - y * other.x
        };
    }

    // Produit scalaire
    double dot(const Vector3D& other) const {
        return x * other.x + y * other.y + z * other.z;
    }
};

#endif