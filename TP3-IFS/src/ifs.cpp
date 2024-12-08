#include "ifs.h"
#include <GL/glut.h>
#include <iostream>


using namespace arma;
using namespace std;

Ifs::Ifs(int level)
{
    m_level = level;
    m_primitive = mat{{-0.5, -0.5, 1.0},
                     { 0.5, -0.5, 1.0},
                     { 0.5,  0.5, 1.0},
                     {-0.5,  0.5, 1.0}};
    
    mat T1 = {{0.5, 0.0, 0.0},
              {0.0, 0.5, 0.0},
              {0.0, 0.0, 1.0}};
    mat T2 = {{1.5, 0.0, 1.5},
              {0.0, 1.5, 1.5},
              {0.0, 0.0, 1.0}};
    
    m_ifs.push_back(T1);
    m_ifs.push_back(T2);

    m_approximation.push_back(m_primitive); 
}

Ifs::~Ifs(void)
{

}

void Ifs::display()
{
    if (m_level < 0 || (size_t)m_level >= m_approximation.size()) {
        return;
    }

    mat points = m_approximation[m_level];

    glColor3f(1.0,1.0,1.0);
    glBegin(GL_POINTS);
    for (uword i = 0; i < points.n_rows; i++) {
        float x = points(i,0);
        float y = points(i,1);
        glVertex2f(x,y);
    }
    glEnd();
}

void Ifs::ComputeApproximation() // il faut peut être mettre des parametres
{
    int level = 20;

    for (int i = 1; i <= level; i++) {
        mat prev = m_approximation[i-1];
        mat current; 
        // On applique toutes les transformations de l'IFS
        for (size_t k = 0; k < m_ifs.size(); k++) {
            mat transformed = prev * m_ifs[k]; 
            // On empile les résultats
            if (current.n_rows == 0)
                current = transformed;
            else
                current = join_cols(current, transformed);
        }
        m_approximation.push_back(current);
    }
}
