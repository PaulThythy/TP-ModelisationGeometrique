#include <iostream>
#include <stdlib.h>
#include <GL/glut.h>
#include <vector>
#include <sstream>
#include <armadillo>

#include <string>

using namespace std;
void affichage(void);

void clavier(unsigned char touche, int x, int y);
void affiche_repere(void);

void mouse(int, int, int, int);
void mouseMotion(int, int);
// void reshape(int,int);
float t = .5;

vector<vector<double>> curve = {};
double steps = 20;
double tSteps = 1/steps;

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
  //  glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE|GLUT_RGB);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45.0f, (GLfloat)200 / (GLfloat)200, 0.1f, 10.0f);
  glMatrixMode(GL_MODELVIEW);
  gluLookAt(0., 0., 4., 0., 0., 0., 0., 1., 0.);
}

//------------------------------------------------------

void displayCourbe(void)
{
  arma::mat P(3,3);
  P(0,0) = 0.0; P(0,1) = 0.0; P(0,2) = 0.0;
  P(1,0) = 1.0; P(1,1) = 1.0; P(1,2) = 0.0;
  P(2,0) = 0.0; P(2,1) = 1.0; P(2,2) = 1.0; 

  arma::mat Mb(3,3);
  Mb(0,0) = 1.0; Mb(0,1) = -2.0; Mb(0,2) = 1.0;
  Mb(1,0) = -2.0; Mb(1,1) = 2.0; Mb(1,2) = 0.0;
  Mb(2,0) = 1.0; Mb(2,1) = 0.0; Mb(2,2) = 0.0;

  arma::mat T(1,3);

  double t = 0;

  while(t<=1) {
    double tt = t*t;
    T(0,0) = tt; T(0,1) = t; T(0,2) = 1; 
    arma::mat vArm = T*Mb*P;
    vector<double> v = { vArm(0,0), vArm(0,1), vArm(0,2) };
    curve.push_back(v);
  }

  for(int i = 0; i <= curve.size(); i++) {
    double v[3] = curve[i];
    double v_next[3] = curve[i+1];
    glColor3f(0.58, 0.0, 0.82);
    glBegin(GL_LINES);
      glVertex3d(v[0], v[1], v[2]);
      glVertex3d(v_next[0], v_next[1], v_next[2]);
    glEnd();
  }


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
