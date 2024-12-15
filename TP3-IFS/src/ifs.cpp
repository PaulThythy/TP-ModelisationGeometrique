#include "ifs.h"
#include <GL/glut.h>
#include <iostream>


using namespace arma;
using namespace std;

Ifs::Ifs()
{
    //TRIANGLE DE SIERPINSKY
    std::vector<arma::mat> transforms = {
        arma::mat{{1.0, .5, .5},
                  {0.0, 0.5, 0.0},
                  {0, 0.0, 0.5}},
        arma::mat{{0.5, 0.0, 0.0},
                  {0.5, 1.0, 0.5},
                  {0.0, 0.0, 0.5}},
        arma::mat{{0.5, 0.0, 0.0},
                  {0.0, 0.5, 0.0},
                  {0.5, 0.5, 1.0}}
    };
    m_transforms = transforms;

    m_primitive = {
        {1, 0, 0}, 
        {0, 1, 0}, 
        {0, 0, 1}
    };

    m_controlPoints = {
        {-1, 0, 1},
        {0, 1, 0},
        {0, 0, 0}
    };

    //ENSEMBLE DE CANTOR
    /*std::vector<arma::mat> transforms = {
        arma::mat{{1.0, 0.6666},
                  {0.0, 0.3333}},
        arma::mat{{0.3333, 0.0},
                  {0.6666, 1.0}}
    };
    m_transforms = transforms;

    m_primitive = {
        {1, 0}, 
        {0, 1}
    };

    m_controlPoints = {
        {0.0, 1.0},
        {0.0, 1.0}
    };*/

    //COURBE DE BEZIER
    /*std::vector<arma::mat> transforms = {
        arma::mat{{1.0, 0.5, 0.25},
                  {0.0, 0.5, 0.5},
                  {0.0, 0.0, 0.25}},
        arma::mat{{0.25, 0.0, 0.0},
                  {0.5, 0.5, 0.0},
                  {0.25, 0.5, 1.0}}
    };
    m_transforms = transforms;

    m_primitive = {
        {1, 0, 0}, 
        {0, 1, 0},
        {0, 0, 1}
    };

    m_controlPoints = {
        {-1, 0, 1},
        {0, 1, 0},
        {0, 0, 0}
    };*/

    //PYRAMIDE DE SIERPINSKY
    /*m_primitive = {
        {
            {0, 1, 0.5},
            {0, 0, 0},
            {0, 0, 0.5}
        },
        {
            {0, 0.5, 0.5},
            {0, 0, 1},
            {0, 0.5, 0}
        },
        {
            {1, 0.5, 0.5},
            {0, 1, 0},
            {0, 0, 0.5}
        },
        {
            {0, 0.5, 1},
            {0, 1, 0},
            {0, 0, 0}
        }
    };*/

    /*std::vector<arma::mat> transforms = {
        arma::mat{{1.0, 0.5, 0.5, 0.5},   // T1
              {0.0, 0.5, 0.0, 0.0},
              {0.0, 0.0, 0.5, 0.0},
              {0.0, 0.0, 0.0, 0.5}},

        arma::mat{{0.5, 0.0, 0.0, 0.0},   // T2
                {0.5, 1.0, 0.5, 0.5},
                {0.0, 0.0, 0.5, 0.0},
                {0.0, 0.0, 0.0, 0.5}},

        arma::mat{{0.5, 0.0, 0.0, 0.0},   // T3
                {0.0, 0.5, 0.0, 0.0},
                {0.5, 0.5, 1.0, 0.5},
                {0.0, 0.0, 0.0, 0.5}},

        arma::mat{{0.5, 0.0, 0.0, 0.0},   // T4
                {0.0, 0.5, 0.0, 0.0},
                {0.0, 0.0, 0.5, 0.0},
                {0.5, 0.5, 0.5, 1.0}}
    };
    m_transforms = transforms;

    m_primitive = {
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    };

    m_controlPoints = {
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    };*/
}

Ifs::~Ifs(void)
{

}

void Ifs::display()
{
    //TRIANGLE DE SIERPINSKY
    glColor3f(1, 1, 1);
    glBegin(GL_TRIANGLES);

    for(int i = 0; i < m_approximation.size(); ++i)
    {
        arma::colvec colonne = m_approximation[i].col(0);
        glVertex3f(colonne(0), colonne(1),colonne(2));
        colonne = m_approximation[i].col(1);
        glVertex3f(colonne(0), colonne(1),colonne(2));
        colonne = m_approximation[i].col(2);
        glVertex3f(colonne(0), colonne(1),colonne(2));
    }
    glEnd();

    //ENSEMBLE DE CANTOR
    /*glColor3f(1, 1, 1);
    glBegin(GL_LINES);

    for(int i = 0; i < m_approximation.size(); ++i)
    {
        arma::colvec colonne = m_approximation[i].col(0);
        glVertex2f(colonne(0), colonne(1));
        colonne = m_approximation[i].col(1);
        glVertex2f(colonne(0), colonne(1));
    }
    glEnd();*/

    //COURBE DE BEZIER
    /*glColor3f(1, 1, 1);
    glBegin(GL_LINES);

    for(int i = 0; i < m_approximation.size(); ++i)
    {
        arma::colvec colonne = m_approximation[i].col(0);
        glVertex2f(colonne(0), colonne(1));
        colonne = m_approximation[i].col(1);
        glVertex2f(colonne(0), colonne(1));
    }
    glEnd();*/

    //PYRAMIDE DE SIERPINSKY
    /*glColor3f(1, 1, 1);
    glBegin(GL_TRIANGLES);

    for (size_t i = 0; i < m_approximation.size(); ++i)
    {
        arma::colvec v0 = m_approximation[i].col(0);
        arma::colvec v1 = m_approximation[i].col(1);
        arma::colvec v2 = m_approximation[i].col(2);
        arma::colvec v3 = m_approximation[i].col(3);

        // Face 1 : (v0, v1, v2)
        glVertex3f((GLfloat)v0(0), (GLfloat)v0(1), (GLfloat)v0(2));
        glVertex3f((GLfloat)v1(0), (GLfloat)v1(1), (GLfloat)v1(2));
        glVertex3f((GLfloat)v2(0), (GLfloat)v2(1), (GLfloat)v2(2));

        // Face 2 : (v0, v1, v3)
        glVertex3f((GLfloat)v0(0), (GLfloat)v0(1), (GLfloat)v0(2));
        glVertex3f((GLfloat)v1(0), (GLfloat)v1(1), (GLfloat)v1(2));
        glVertex3f((GLfloat)v3(0), (GLfloat)v3(1), (GLfloat)v3(2));

        // Face 3 : (v0, v2, v3)
        glVertex3f((GLfloat)v0(0), (GLfloat)v0(1), (GLfloat)v0(2));
        glVertex3f((GLfloat)v2(0), (GLfloat)v2(1), (GLfloat)v2(2));
        glVertex3f((GLfloat)v3(0), (GLfloat)v3(1), (GLfloat)v3(2));

        // Face 4 : (v1, v2, v3)
        glVertex3f((GLfloat)v1(0), (GLfloat)v1(1), (GLfloat)v1(2));
        glVertex3f((GLfloat)v2(0), (GLfloat)v2(1), (GLfloat)v2(2));
        glVertex3f((GLfloat)v3(0), (GLfloat)v3(1), (GLfloat)v3(2));
    }

    glEnd();*/
}

void Ifs::ComputeApproximation(int iteration)
{
    m_approximation.resize(0);
    m_approximation.push_back(m_primitive);

    int nbTransforms = m_transforms.size();

    //chaque itération ajoute une puissance de 2 de primitives
    size_t nbMaxTransforms = m_approximation.size();
    for (int i = 0; i < iteration; ++i) {
        nbMaxTransforms *= m_transforms.size();
    }

    int n = sqrt(m_primitive.size());

    m_approximation.clear();
    m_approximation.resize(nbMaxTransforms);
    //matrice identité à l'itération 0
    m_approximation[0] = eye<arma::mat>(n, n);

    int size = 1;

    for (int i = 0; i < iteration; ++i)
    {
        //Pour chaque transform
        //boucle à l'envers car : optimisation pour pas utiliser trop de mémoire
        //rempli de gauche à droite, mais par itération ça rempli de droite à gauche
        for (int j = nbTransforms-1; j >= 0; j--)
        {
            //Chaque feuille sera multiplié par la transform j
            //size le nombre de primitives 
            for (int k = 0; k < size; k++)
            {
                mat approx = m_transforms[j] * m_approximation[k];
                m_approximation[j * size + k] = approx;
            }
        }
        //MAJ du nombre de feuilles
        size = size * nbTransforms;
    }

    //espace barycentrique vers l'espace modélisation
    for (int i = 0; i < size; ++i)
    {
        m_approximation[i] = m_controlPoints * m_approximation[i];
    }
}

void Ifs::printApproximation() 
{
    for (size_t i = 0; i < m_approximation.size(); ++i) {
        std::cout << "Matrice " << i + 1 << ":" << std::endl;
        m_approximation[i].print();
        std::cout << std::endl;
    }
}
