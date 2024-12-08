#include <iostream>
#include <stdlib.h>
#include <GL/glut.h>
#include <vector>
#include <sstream>
#include <string>
#include <cmath>
#include <armadillo>

using namespace std;
void affichage(void);

void clavier(unsigned char touche, int x, int y);
void affiche_repere(void);

void mouse(int, int, int, int);
void mouseMotion(int, int);
// void reshape(int,int);
float t = .5;

// Structure pour représenter le repère de Frenet
struct FrenetFrame
{
  arma::mat T; // Vecteur tangent
  arma::mat N; // Vecteur normal
  arma::mat B; // Vecteur binormal
};

// variables globales pour OpenGL
bool mouseLeftDown;
bool mouseRightDown;
bool mouseMiddleDown;
float mouseX, mouseY;
float cameraAngleX;
float cameraAngleY;
float cameraDistance = 0.;

const int width = 600;
const int height = 600;

// Degrés des NURBS en u et v
int p = 3; // Degré en u
int q = 3; // Degré en v

// Nombre de points de contrôle en u et v
int n = 3; // n+1 = 4 points de contrôle en u
int m = 3; // m+1 = 4 points de contrôle en v

GLfloat ctrlPoints[4][4][3] = {
    {{-1.0, -1.0, 3.0}, {-0.5, -1.0, 2.0}, {0.5, -1.0, -1.0}, {1.0, -1.0, 2.0}},
    {{-1.0, -0.5, 1.0}, {-0.5, -0.5, 3.0}, {0.5, -0.5, 0.0}, {1.0, -0.5, -1.0}},
    {{-1.0, 0.5, 3.0}, {-0.5, 0.5, 0.0}, {0.5, 0.5, 3.0}, {1.0, 0.5, 3.0}},
    {{-1.0, 1.0, -2.0}, {-0.5, 1.0, -2.0}, {0.5, 1.0, 0.0}, {1.0, 1.0, -1.0}}};

GLfloat knotU[] = {0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 3.0, 3.0};
GLfloat knotV[] = {0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 3.0, 3.0};

GLfloat u = (knotU[p] + knotU[n+1]) * 0.5f;
GLfloat v = (knotV[q] + knotV[m+1]) * 0.5f; 

GLfloat weights[4][4] = {
    {1.0, 1.0, 1.0, 1.0},
    {1.0, 1.0, 1.0, 1.0},
    {1.0, 1.0, 1.0, 1.0},
    {1.0, 1.0, 1.0, 1.0}};

// Trouver l'indice du noeud (knot span) pour un paramètre donné
int findKnotSpan(int n, int p, float u, GLfloat knotVector[])
{
  if (u >= knotVector[n + 1])
  {
    return n;
  }
  if (u <= knotVector[p])
  {
    return p;
  }
  int low = p;
  int high = n + 1;
  int mid = (low + high) / 2;
  while (u < knotVector[mid] || u >= knotVector[mid + 1])
  {
    if (u < knotVector[mid])
    {
      high = mid;
    }
    else
    {
      low = mid;
    }
    mid = (low + high) / 2;
  }
  return mid;
}

// Calculer les fonctions de base B-spline non nulles
void basisFunctions(int i, float u, int p, GLfloat knotVector[], float N[])
{
  N[0] = 1.0;
  float left[p + 1];
  float right[p + 1];
  for (int j = 1; j <= p; j++)
  {
    left[j] = u - knotVector[i + 1 - j];
    right[j] = knotVector[i + j] - u;
    float saved = 0.0;
    for (int r = 0; r < j; r++)
    {
      float temp = N[r] / (right[r + 1] + left[j - r]);
      N[r] = saved + right[r + 1] * temp;
      saved = left[j - r] * temp;
    }
    N[j] = saved;
  }
}

// Évaluer le point de surface NURBS à (u, v)
arma::mat surfacePoint(float u, float v)
{
  int numCtrlPointsU = n + 1; // n = degré en u
  int numCtrlPointsV = m + 1; // m = degré en v

  // Trouver les indices des nœuds pour u et v
  int spanU = findKnotSpan(n, p, u, knotU);
  int spanV = findKnotSpan(m, q, v, knotV);

  // Calculer les fonctions de base pour u et v
  float Nu[p + 1];
  float Nv[q + 1];
  basisFunctions(spanU, u, p, knotU, Nu);
  basisFunctions(spanV, v, q, knotV, Nv);

  arma::mat S = arma::zeros<arma::mat>(3, 1);
  float denominator = 0.0f;

  // Calculer le point de surface en utilisant les fonctions de base et les poids
  for (int l = 0; l <= q; l++)
  {
    float temp = 0.0f;
    arma::mat tempPoint = arma::zeros<arma::mat>(3, 1);
    for (int k = 0; k <= p; k++)
    {
      int uIndex = spanU - p + k;
      int vIndex = spanV - q + l;
      float w = weights[uIndex][vIndex];
      float basis = Nu[k] * Nv[l] * w;
      arma::mat P(3, 1);
      P(0, 0) = ctrlPoints[uIndex][vIndex][0];
      P(1, 0) = ctrlPoints[uIndex][vIndex][1];
      P(2, 0) = ctrlPoints[uIndex][vIndex][2];
      tempPoint += basis * P;
      denominator += basis;
    }
    S += tempPoint;
  }
  S /= denominator;
  return S;
}

arma::mat computePartialDerivativeU(float u, float v)
{
  float delta = 0.01f; // Pas de calcul
  arma::mat S1 = surfacePoint(u + delta, v);
  arma::mat S0 = surfacePoint(u - delta, v);
  return (S1 - S0) / (2.0f * delta);
}

arma::mat computePartialDerivativeV(float u, float v)
{
  float delta = 0.01f; // Pas de calcul
  arma::mat S1 = surfacePoint(u, v + delta);
  arma::mat S0 = surfacePoint(u, v - delta);
  return (S1 - S0) / (2.0f * delta);
}

arma::mat computeSecondPartialDerivativeUU(float u, float v)
{
  float delta = 0.01f; // Pas de calcul
  arma::mat S1 = surfacePoint(u + delta, v);
  arma::mat S0 = surfacePoint(u, v);
  arma::mat S2 = surfacePoint(u - delta, v);
  return (S1 - 2.0f * S0 + S2) / (delta * delta);
}

arma::mat computeSecondPartialDerivativeVV(float u, float v)
{
  float delta = 0.01f; // Pas de calcul
  arma::mat S1 = surfacePoint(u, v + delta);
  arma::mat S0 = surfacePoint(u, v);
  arma::mat S2 = surfacePoint(u, v - delta);
  return (S1 - 2.0f * S0 + S2) / (delta * delta);
}

float computeCurvatureRadiusAt(float u, float v)
{
  arma::mat du = computePartialDerivativeU(u, v);
  arma::mat dv = computePartialDerivativeV(u, v);
  arma::mat duu = computeSecondPartialDerivativeUU(u, v);
  arma::mat dvv = computeSecondPartialDerivativeVV(u, v);

  float numerator = arma::norm(arma::cross(duu, dvv));
  float denominator = std::pow(arma::norm(du), 3);
  return numerator / denominator;
}

FrenetFrame computeFrenetFrameAt(float u, float v)
{
  // Dérivées premières
  arma::mat du = computePartialDerivativeU(u, v);
  arma::mat dv = computePartialDerivativeV(u, v);

  // Tangente (vecteur T)
  arma::mat T = arma::normalise(du + dv);

  // Dérivées secondes
  arma::mat duu = computeSecondPartialDerivativeUU(u, v);
  arma::mat dvv = computeSecondPartialDerivativeVV(u, v);

  // Normale (vecteur N)
  arma::mat N = arma::normalise(duu + dvv);

  // Binormale (vecteur B)
  arma::mat B = arma::cross(T, N);

  return {T, N, B};
}

void displayFrenetFrameAt(float u, float v)
{
  FrenetFrame frame = computeFrenetFrameAt(u, v);
  arma::mat point = surfacePoint(u, v);

  GLfloat pointArray[3] = {
      static_cast<GLfloat>(point(0, 0)),
      static_cast<GLfloat>(point(1, 0)),
      static_cast<GLfloat>(point(2, 0))};

  glBegin(GL_LINES);

  // Tangente T
  arma::mat tangent = point + frame.T; // Evaluation explicite
  GLfloat tangentArray[3] = {
      static_cast<GLfloat>(tangent(0, 0)),
      static_cast<GLfloat>(tangent(1, 0)),
      static_cast<GLfloat>(tangent(2, 0))};
  glColor3f(1.0f, 0.0f, 0.0f); // Rouge
  glVertex3fv(pointArray);
  glVertex3fv(tangentArray);

  // Normale N
  arma::mat normal = point + frame.N; // Evaluation explicite
  GLfloat normalArray[3] = {
      static_cast<GLfloat>(normal(0, 0)),
      static_cast<GLfloat>(normal(1, 0)),
      static_cast<GLfloat>(normal(2, 0))};
  glColor3f(0.0f, 1.0f, 0.0f); // Vert
  glVertex3fv(pointArray);
  glVertex3fv(normalArray);

  // Binormale B
  arma::mat binormal = point + frame.B; // Evaluation explicite
  GLfloat binormalArray[3] = {
      static_cast<GLfloat>(binormal(0, 0)),
      static_cast<GLfloat>(binormal(1, 0)),
      static_cast<GLfloat>(binormal(2, 0))};
  glColor3f(0.0f, 0.0f, 1.0f); // Bleu
  glVertex3fv(pointArray);
  glVertex3fv(binormalArray);

  glEnd();
}

void initOpenGl()
{

  // lumiere
  glClearColor(.5, .5, 0.5, 0.0);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  GLfloat l_pos[] = {
      3.,
      3.5,
      3.0,
      1.0};
  glLightfv(GL_LIGHT0, GL_POSITION, l_pos);

  glLightfv(GL_LIGHT0, GL_DIFFUSE, l_pos);
  glLightfv(GL_LIGHT0, GL_SPECULAR, l_pos);
  glEnable(GL_COLOR_MATERIAL);

  glDepthFunc(GL_LESS);
  glEnable(GL_DEPTH_TEST);
  // glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
  //  glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE|GLUT_RGB);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45.0f, (GLfloat)width / (GLfloat)height, 0.1f, 10000.0f);
  glMatrixMode(GL_MODELVIEW);
  gluLookAt(0., 0., 4., 0., 0., 0., 0., 1., 0.);
}

//------------------------------------------------------

void displaySurface(void)
{
  GLUnurbsObj *nurbs = gluNewNurbsRenderer();

  gluNurbsProperty(nurbs, GLU_DISPLAY_MODE, GLU_OUTLINE_POLYGON);

  gluBeginSurface(nurbs);
  gluNurbsSurface(nurbs,
                  8, knotU,
                  8, knotV,
                  4 * 3,
                  3,
                  &ctrlPoints[0][0][0],
                  4,
                  4,
                  GL_MAP2_VERTEX_3);
  gluEndSurface(nurbs);

  gluDeleteNurbsRenderer(nurbs);
}

int main(int argc, char **argv)
{
  /* initialisation de glut et creation de la fenetre */
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGB);
  glutInitWindowPosition(200, 200);
  glutInitWindowSize(width, height);
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
  displaySurface();
  displayFrenetFrameAt(u, v);
  glPopMatrix();
  glFlush();
  glutSwapBuffers();
}

//------------------------------------------------------

//------------------------------------------------------
void clavier(unsigned char touche, int x, int y)
{

  switch (touche)
  {
  case '+': //
    t += .1;
    if (t > 1)
      t = 1;
    glutPostRedisplay();
    break;
  case '-': //* ajustement du t
    t -= .1;
    if (t < 0)
      t = 0;
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
  case 'u':
    u += 0.05f;
    if (u > knotU[n+1])
        u = knotU[n+1];
    if (u < knotU[p])
        u = knotU[p];
    glutPostRedisplay();
    break;
  case 'U':
    u -= 0.05f;
    if (u > knotU[n+1])
        u = knotU[n+1];
    if (u < knotU[p])
        u = knotU[p];
    glutPostRedisplay();
    break;
  case 'v':
    v += 0.05f;
    if (v > knotV[m+1])
        v = knotV[m+1];
    if (v < knotV[q])
        v = knotV[q];
    glutPostRedisplay();
    break;
  case 'V':
    v -= 0.05f;
    if (v > knotV[m+1])
        v = knotV[m+1];
    if (v < knotV[q])
        v = knotV[q];
    glutPostRedisplay();
    break;
  case 'q':
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