#ifndef __IFS
#define __IFS

#define NB_TRANSFOS 3
#include <armadillo>

#define BARY 1
#define HOMOGENE 2
#include <GL/glut.h>

using namespace arma;
using namespace std;

class Ifs
{
public:
	vector<mat> mIfs;
	mat mPrimitive;
	mat mControlPoints;
	vector<mat> mApproximation;

	Ifs(void);
	~Ifs(void);
	void display(int level);
	void ComputeApproximation(); // il faut peut être mettre des parametres
};

#endif
