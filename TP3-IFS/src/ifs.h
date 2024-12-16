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
	vector<mat> m_transforms;
	mat m_primitive;
	mat m_controlPoints;
	vector<mat> m_approximation;

	Ifs(void);
	~Ifs(void);
	void display();
	void ComputeApproximation(int iteration);
};

#endif
