#include "ifs.h"
#include <GL/glut.h>
#include <iostream>
// #include "../armadillo/include/armadillo"


using namespace arma;
using namespace std;

Ifs::Ifs(void)
{
    //triangle en coordonnées barycentriques
    mPrimitive = {
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}
    };

    //points de contrôles en coordonnées géométriques
    mControlPoints = {
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}
    };

    mIfs = {
        
    };
}

Ifs::~Ifs(void)
{

}

void Ifs::display(int level)
{

}

void Ifs::ComputeApproximation() // il faut peut être mettre des parametres
{

}
