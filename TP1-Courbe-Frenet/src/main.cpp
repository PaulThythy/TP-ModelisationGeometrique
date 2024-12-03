#include <iostream>
#include <stdlib.h>
#include <GL/glut.h>
#include <vector>
#include <sstream>
#include <armadillo>
#include <string>

#include "Vector3D.h"

using namespace std;
void affichage(void);

void clavier(unsigned char touche, int x, int y);
void affiche_repere(void);

void mouse(int, int, int, int);
void mouseMotion(int, int);
// void reshape(int,int);

// Structure pour représenter le repère de Frenet
struct FrenetFrame
{
  Vector3D T; // Vecteur tangent
  Vector3D N; // Vecteur normal
  Vector3D B; // Vecteur binormal
};

// Points de contrôle, degré et vecteur nodal
std::vector<Vector3D> pointsControle = {
    {-1.5, -1.0, 0.0},
    {-0.5, 1.0, 0.0},
    {0.5, -1.0, 0.0},
    {1.5, 1.0, 0.0},
    {2.5, -1.0, 0.0}};
int degre = 3;
arma::vec vecteurNodal = {0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 4};

// Stockage des points de la courbe
std::vector<Vector3D> pointsCourbe;

int resolution = 10;
// Incrément
double deltaU = 1.0 / resolution; // Par exemple, si résolution = 1000
double current_u = 0.5f;

// Variables globales pour u_min et u_max
double u_min;
double u_max;

// variables globales pour OpenGL
bool mouseLeftDown;
bool mouseRightDown;
bool mouseMiddleDown;
float mouseX, mouseY;
float cameraAngleX;
float cameraAngleY;
float cameraDistance = 0.;

// constantes pour les materieux
float no_mat[] = {0.0f, 0.0f, 0.0f, 1.0f};
float mat_ambient[] = {0.7f, 0.7f, 0.7f, 1.0f};
float mat_ambient_color[] = {0.8f, 0.8f, 0.2f, 1.0f};
float mat_diffuse[] = {0.1f, 0.5f, 0.8f, 1.0f};
float mat_specular[] = {1.0f, 1.0f, 1.0f, 1.0f};
float no_shininess = 0.0f;
float low_shininess = 5.0f;
float high_shininess = 100.0f;
float mat_emission[] = {0.3f, 0.2f, 0.2f, 0.0f};

// Fonction de base B-Spline récursive
double B_Spline(int i, int d, double u, const arma::vec &knots)
{
  if (d == 0)
  {
    if (knots[i] <= u && u < knots[i + 1])
      return 1.0;
    else
      return 0.0;
  }
  else
  {
    double left = 0.0, right = 0.0;
    double denom1 = knots[i + d] - knots[i];
    if (denom1 != 0)
      left = (u - knots[i]) / denom1 * B_Spline(i, d - 1, u, knots);
    double denom2 = knots[i + d + 1] - knots[i + 1];
    if (denom2 != 0)
      right = (knots[i + d + 1] - u) / denom2 * B_Spline(i + 1, d - 1, u, knots);
    return left + right;
  }
}

// Fonction pour générer les points de la courbe B-Spline
std::vector<Vector3D> genererBSpline(const std::vector<Vector3D> &controlePoints, int d, const arma::vec &knots, int resolution)
{
  std::vector<Vector3D> courbe;
  int n = controlePoints.size() - 1;
  double u_min = knots[d];
  double u_max = knots[n + 1];

  for (int i = 0; i <= resolution; ++i)
  {
    double u = u_min + (u_max - u_min) * i / resolution;
    double x = 0.0, y = 0.0, z = 0.0;
    for (int j = 0; j <= n; ++j)
    {
      double Nj = B_Spline(j, d, u, knots);
      x += Nj * controlePoints[j].x;
      y += Nj * controlePoints[j].y;
      z += Nj * controlePoints[j].z;
    }
    courbe.push_back(Vector3D{x, y, z});
  }
  return courbe;
}

// Fonction pour calculer le point sur la courbe à un paramètre u donné
Vector3D computePointAt(double u)
{
  Vector3D point = {0.0, 0.0, 0.0};
  int n = pointsControle.size() - 1;
  for (int i = 0; i <= n; ++i)
  {
    double Ni = B_Spline(i, degre, u, vecteurNodal);
    point = point + pointsControle[i] * Ni;
  }
  return point;
}

// Fonction pour calculer la dérivée première par différences finies
Vector3D computeFirstDerivativeAt(double u)
{
  double du = 1e-5; // petit incrément
  double u_forward = std::min(u + du, u_max);
  double u_backward = std::max(u - du, u_min);
  Vector3D point_forward = computePointAt(u_forward);
  Vector3D point_backward = computePointAt(u_backward);
  return (point_forward - point_backward) / (2.0 * du);
}

// Fonction pour calculer la dérivée seconde par différences finies
Vector3D computeSecondDerivativeAt(double u)
{
  double du = 1e-5; // petit incrément
  double u_forward = std::min(u + du, u_max);
  double u_backward = std::max(u - du, u_min);
  Vector3D point_forward = computePointAt(u_forward);
  Vector3D point_current = computePointAt(u);
  Vector3D point_backward = computePointAt(u_backward);
  return (point_forward - point_current * 2.0 + point_backward) / (du * du);
}

// Fonction pour calculer le repère de Frenet à un paramètre u donné
FrenetFrame computeFrenetFrameAt(double u)
{
  Vector3D dC = computeFirstDerivativeAt(u);
  Vector3D ddC = computeSecondDerivativeAt(u);
  Vector3D T = dC.normalize();
  Vector3D numerator = ddC - T * (ddC.dot(T));
  Vector3D N = numerator.normalize();
  Vector3D B = T.cross(N);
  return FrenetFrame{ T, N, B };
}

// Fonction pour calculer le rayon de courbure à un paramètre u donné
double computeCurvatureAt(double u)
{
  Vector3D dC = computeFirstDerivativeAt(u);
  Vector3D ddC = computeSecondDerivativeAt(u);
  Vector3D crossProduct = dC.cross(ddC);
  double numerator = crossProduct.norm();
  double denominator = pow(dC.norm(), 3);
  double curvature = numerator / denominator;
  double radius = (curvature != 0.0) ? 1.0 / curvature : 0.0;
  return radius;
}

// Initialisation des points de la courbe
void initialiserCourbe()
{
  int n = pointsControle.size() - 1;
  u_min = vecteurNodal[degre];
  u_max = vecteurNodal[n + 1];
  pointsCourbe = genererBSpline(pointsControle, degre, vecteurNodal, resolution);
  current_u = u_min; // Initialiser current_u à u_min
}

// Affiche un cercle le long de la courbe
void displayCircle() {
  int nCircleSegments = 20; // Nombre de segments du cercle
    double circleRadius = 0.2; // Rayon du cercle
    const double EPSILON = 1e-6;
    double deltaU_circle = (u_max - u_min) / (resolution * 10); // Ajustez le facteur selon vos besoins

    glColor3f(0.7f, 0.7f, 0.7f); // Couleur des cercles

    for (double u = u_min; u <= u_max + EPSILON; u += deltaU_circle)
    {
        Vector3D point = computePointAt(u);
        FrenetFrame frame = computeFrenetFrameAt(u);

        glBegin(GL_LINE_LOOP);
        for (int j = 0; j < nCircleSegments; ++j)
        {
            double theta = 2.0 * M_PI * double(j) / double(nCircleSegments);
            Vector3D circlePoint = point + (frame.N * cos(theta) + frame.B * sin(theta)) * circleRadius;
            glVertex3f(circlePoint.x, circlePoint.y, circlePoint.z);
        }
        glEnd();
    }
}

void initOpenGl()
{

  // lumiere

  glClearColor(.5, .5, 0.5, 0.0);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  GLfloat l_pos[] = {3., 3.5, 3.0, 1.0};
  glLightfv(GL_LIGHT0, GL_POSITION, l_pos);

  glLightfv(GL_LIGHT0, GL_DIFFUSE, l_pos);
  glLightfv(GL_LIGHT0, GL_SPECULAR, l_pos);
  glEnable(GL_COLOR_MATERIAL);

  glDepthFunc(GL_LESS);
  glEnable(GL_DEPTH_TEST);
  // glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
  // glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE|GLUT_RGB);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45.0f, (GLfloat)200 / (GLfloat)200, 0.1f, 10.0f);
  glMatrixMode(GL_MODELVIEW);
  gluLookAt(0., 0., 4., 0., 0., 0., 0., 1., 0.);
}

//------------------------------------------------------

void displayCourbe(void)
{
  if (pointsCourbe.empty())
    return;

  glColor3f(1.0f, 1.0f, 1.0f); // Couleur de la courbe
  glBegin(GL_LINE_STRIP);
  for (const auto &point : pointsCourbe)
  {
    glVertex3f(point.x, point.y, point.z);
  }
  glEnd();

  // Optionnel : Afficher les points de contrôle
  glPointSize(5.0f);
  glColor3f(1.0f, 0.0f, 0.0f); // Couleur des points de contrôle
  glBegin(GL_POINTS);
  for (const auto &point : pointsControle)
  {
    glVertex3f(point.x, point.y, point.z);
  }
  glEnd();
}

int main(int argc, char **argv)
{
  /* initialisation de glut et creation
     de la fenetre */
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGB);
  glutInitWindowPosition(200, 200);
  glutInitWindowSize(600, 600);
  glutCreateWindow("ifs");

  /* Initialisation d'OpenGL */
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glColor3f(1.0, 1.0, 1.0);
  glPointSize(1.0);

  // ifs = new Ifs();
  /* enregistrement des fonctions de rappel */
  glutDisplayFunc(affichage);
  glutKeyboardFunc(clavier);
  glutMouseFunc(mouse);
  glutMotionFunc(mouseMotion);
  //-------------------------------

  //-------------------------------
  initOpenGl();
  //-------------------------------

  initialiserCourbe();

  /* Entree dans la boucle principale glut */
  glutMainLoop();
  return 0;
}
//------------------------------------------------------
void affiche_repere(void)
{
  glBegin(GL_LINES);
  glColor3f(1.0, 0.0, 0.0);
  glVertex2f(0., 0.);
  glVertex2f(1., 0.);
  glEnd();

  glBegin(GL_LINES);
  glColor3f(0.0, 1.0, 0.0);
  glVertex2f(0., 0.);
  glVertex2f(0., 1.);
  glEnd();

  glBegin(GL_LINES);
  glColor3f(0.0, 0.0, 1.0);
  glVertex3f(0., 0., 0.);
  glVertex3f(0., 0., 1.);
  glEnd();
}

//-----------------------------------------------------

//------------------------------------------------------
void affichage(void)
{
  glMatrixMode(GL_MODELVIEW);
  /* effacement de l'image avec la couleur de fond */
  //	glClear(GL_COLOR_BUFFER_BIT);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  //       glClearDepth(10.0f);                         // 0 is near, >0 is far

  glPushMatrix();
  glTranslatef(0, 0, cameraDistance);
  glRotatef(cameraAngleX, 1., 0., 0.);
  glRotatef(cameraAngleY, 0., 1., 0.);
  affiche_repere();
  displayCourbe();
  displayCircle();

  // Calculer et dessiner le repère de Frenet
  Vector3D point = computePointAt(current_u);
  FrenetFrame frenet = computeFrenetFrameAt(current_u);

  // Calculer le rayon de courbure
  double radius = computeCurvatureAt(current_u);

  // Vérifier si le rayon est valide
  if (radius > 0 && radius < std::numeric_limits<double>::infinity())
  {
    // Calculer le centre du cercle osculateur
    Vector3D center = point + frenet.N * radius;

    // Dessiner le cercle osculateur
    int num_segments = 100;   // Nombre de segments pour approximer le cercle
    glColor3f(1.0, 1.0, 0.0); // Couleur du cercle osculateur (jaune)
    glBegin(GL_LINE_LOOP);
    for (int i = 0; i < num_segments; ++i)
    {
      double theta = 2.0 * M_PI * double(i) / double(num_segments);
      Vector3D circle_point = center + (frenet.N * cos(theta) + frenet.B * sin(theta)) * radius;
      glVertex3f(circle_point.x, circle_point.y, circle_point.z);
    }
    glEnd();
  }

  glBegin(GL_LINES);
  // Vecteur tangent T (rouge)
  glColor3f(1.0, 0.0, 0.0);
  glVertex3f(point.x, point.y, point.z);
  glVertex3f(point.x + frenet.T.x, point.y + frenet.T.y, point.z + frenet.T.z);

  // Vecteur normal N (vert)
  glColor3f(0.0, 1.0, 0.0);
  glVertex3f(point.x, point.y, point.z);
  glVertex3f(point.x + frenet.N.x, point.y + frenet.N.y, point.z + frenet.N.z);

  // Vecteur binormal B (bleu)
  glColor3f(0.0, 0.0, 1.0);
  glVertex3f(point.x, point.y, point.z);
  glVertex3f(point.x + frenet.B.x, point.y + frenet.B.y, point.z + frenet.B.z);
  glEnd();

  glPopMatrix();
  /* on force l'affichage du resultat */

  glFlush();
  glutSwapBuffers();
}

//------------------------------------------------------

//------------------------------------------------------
void clavier(unsigned char touche, int x, int y)
{

  switch (touche)
  {
  case 'a': // Augmenter u
    current_u += deltaU;
    if (current_u > u_max)
      current_u = u_max;
    std::cout << "currentU : " << current_u << std::endl;
    glutPostRedisplay();
    break;
  case 'z': // Diminuer u
    current_u -= deltaU;
    if (current_u < u_min)
      current_u = u_min;
    std::cout << "currentU : " << current_u << std::endl;
    glutPostRedisplay();
    break;
  case 'f': //* affichage en mode fil de fer
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glutPostRedisplay();
    break;
  case 'p': //* affichage du carre plein
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glutPostRedisplay();
    break;
  case 's': //* Affichage en mode sommets seuls
    glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
    glutPostRedisplay();
    break;

  case 'q': //*la touche 'q' permet de quitter le programme
    exit(0);
  }
}
void mouse(int button, int state, int x, int y)
{
  mouseX = x;
  mouseY = y;

  if (button == GLUT_LEFT_BUTTON)
  {
    if (state == GLUT_DOWN)
    {
      mouseLeftDown = true;
    }
    else if (state == GLUT_UP)
      mouseLeftDown = false;
  }

  else if (button == GLUT_RIGHT_BUTTON)
  {
    if (state == GLUT_DOWN)
    {
      mouseRightDown = true;
    }
    else if (state == GLUT_UP)
      mouseRightDown = false;
  }

  else if (button == GLUT_MIDDLE_BUTTON)
  {
    if (state == GLUT_DOWN)
    {
      mouseMiddleDown = true;
    }
    else if (state == GLUT_UP)
      mouseMiddleDown = false;
  }
}

void mouseMotion(int x, int y)
{
  if (mouseLeftDown)
  {
    cameraAngleY += (x - mouseX);
    cameraAngleX += (y - mouseY);
    mouseX = x;
    mouseY = y;
  }
  if (mouseRightDown)
  {
    cameraDistance += (y - mouseY) * 0.2f;
    mouseY = y;
  }

  glutPostRedisplay();
}
