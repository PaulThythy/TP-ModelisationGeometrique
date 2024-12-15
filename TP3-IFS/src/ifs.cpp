#include "ifs.h"
#include <GL/glut.h>
#include <iostream>


using namespace arma;
using namespace std;

Ifs::Ifs(int level)
{
    m_level = level;

    // Triangle équilatéral en coordonnées cartésiennes homogènes
    double sqrt3 = std::sqrt(3.0);
    m_primitive = {
        {0.0,     0.0,     1.0},  // Sommet A
        {1.0,     0.0,     1.0},  // Sommet B
        {0.5, sqrt3 / 2.0, 1.0}   // Sommet C
    };

    // Transformations du triangle de Sierpiński en cartésien homogène
    mat T1 = {
        {0.5, 0.0, 0.0},
        {0.0, 0.5, 0.0},
        {0.0, 0.0, 1.0}
    };

    mat T2 = {
        {0.5, 0.0, 0.5},
        {0.0, 0.5, 0.0},
        {0.0, 0.0, 1.0}
    };

    mat T3 = {
        {0.5,    0.0,   0.25},
        {0.0,    0.5, sqrt3 / 4.0},
        {0.0,    0.0,   1.0}
    };

    // Ajouter les transformations au vecteur m_ifs
    m_ifs.push_back(T1);
    m_ifs.push_back(T2);
    m_ifs.push_back(T3);

    // Initialisation du vecteur d'approximation avec la primitive
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
    std::cout << points << std::endl;

    glColor3f(1.0,0.0,0.0);
    glBegin(GL_POINTS);
    glPointSize(3.0f);
    for (uword i = 0; i < points.n_rows; i++) {
        float x = points(i,0);
        float y = points(i,1);
        glVertex2f(x,y);
    }
    glEnd();
}

void Ifs::ComputeApproximation() // il faut peut être mettre des parametres
{
    for (int i = 1; i <= m_level; i++) {
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
    
    // mise à l'échelle pour rendre les points plus visibles
    for (auto &matPoints : m_approximation) {
        matPoints.col(0) *= 2.0;
        matPoints.col(1) *= 2.0;
    }
}
