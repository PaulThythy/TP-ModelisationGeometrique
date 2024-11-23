#include <iostream>
#include <stdlib.h>
#include <GL/glut.h>
#include <vector>
#include <sstream>
#include <string>
#include <cmath>

using namespace std;
void affichage(void);

void clavier(unsigned char touche, int x, int y);
void affiche_repere(void);

void mouse(int, int, int, int);
void mouseMotion(int, int);
//void reshape(int,int);
float t = .5;

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

GLfloat u = 0.1;
GLfloat v = 0.1;

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
    {{-1.0, 1.0, -2.0}, {-0.5, 1.0, -2.0}, {0.5, 1.0, 0.0}, {1.0, 1.0, -1.0}}
};

GLfloat knotU[] = {0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 3.0, 3.0};
GLfloat knotV[] = {0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 3.0, 3.0};

GLfloat weights[4][4] = {
    {1.0, 1.0, 1.0, 1.0},
    {1.0, 1.0, 1.0, 1.0},
    {1.0, 1.0, 1.0, 1.0},
    {1.0, 1.0, 1.0, 1.0}
};

// constantes pour les materieux
float no_mat[] = {
  0.0f,
  0.0f,
  0.0f,
  1.0f
};
float mat_ambient[] = {
  0.7f,
  0.7f,
  0.7f,
  1.0f
};
float mat_ambient_color[] = {
  0.8f,
  0.8f,
  0.2f,
  1.0f
};
float mat_diffuse[] = {
  0.1f,
  0.5f,
  0.8f,
  1.0f
};
float mat_specular[] = {
  1.0f,
  1.0f,
  1.0f,
  1.0f
};
float no_shininess = 0.0f;
float low_shininess = 5.0f;
float high_shininess = 100.0f;
float mat_emission[] = {
  0.3f,
  0.2f,
  0.2f,
  0.0f
};

void initOpenGl() {

  //lumiere 
  glClearColor(.5, .5, 0.5, 0.0);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  GLfloat l_pos[] = {
    3.,
    3.5,
    3.0,
    1.0
  };
  glLightfv(GL_LIGHT0, GL_POSITION, l_pos);

  glLightfv(GL_LIGHT0, GL_DIFFUSE, l_pos);
  glLightfv(GL_LIGHT0, GL_SPECULAR, l_pos);
  glEnable(GL_COLOR_MATERIAL);

  glDepthFunc(GL_LESS);
  glEnable(GL_DEPTH_TEST);
  //glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
  // glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE|GLUT_RGB);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45.0f, (GLfloat) width / (GLfloat) height, 0.1f, 10000.0f);
  glMatrixMode(GL_MODELVIEW);
  gluLookAt(0., 0., 4., 0., 0., 0., 0., 1., 0.);

}

//------------------------------------------------------

// Fonction pour calculer les fonctions de base B-spline N_{i,p}(u)
double Nip(int i, int p, double u, const GLfloat* U) {
  if (p == 0) {
      if (U[i] <= u && u < U[i+1]) {
          return 1.0;
      } else {
          return 0.0;
      }
  } else {
      double denom1 = U[i+p] - U[i];
      double denom2 = U[i+p+1] - U[i+1];
      double term1 = 0.0, term2 = 0.0;

      if (denom1 != 0.0) {
          term1 = (u - U[i]) / denom1 * Nip(i, p - 1, u, U);
      }
      if (denom2 != 0.0) {
          term2 = (U[i+p+1] - u) / denom2 * Nip(i + 1, p - 1, u, U);
      }
      return term1 + term2;
  }
}

// Fonction pour calculer la dérivée des fonctions de base B-spline N'_{i,p}(u)
double NipDerivative(int i, int p, double u, const GLfloat* U) {
  if (p == 0) {
      return 0.0;
  } else {
      double denom1 = U[i+p] - U[i];
      double denom2 = U[i+p+1] - U[i+1];
      double term1 = 0.0, term2 = 0.0;

      if (denom1 != 0.0) {
          term1 = (1.0 / denom1) * Nip(i, p - 1, u, U);
      }
      if (denom2 != 0.0) {
          term2 = (-1.0 / denom2) * Nip(i + 1, p - 1, u, U);
      }
      return term1 + term2 + 
          (u - U[i]) / denom1 * NipDerivative(i, p - 1, u, U) + 
          (U[i+p+1] - u) / denom2 * NipDerivative(i + 1, p - 1, u, U);
  }
}

void displayCourbe(void) {
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

void afficheVecteursTangents(GLfloat u, GLfloat v) {
  GLfloat point[3], du[3], dv[3], normale[3];

  GLfloat delta = 0.01f;
  GLfloat point_u_plus_delta[3], point_v_plus_delta[3];

  double numerator[3] = {0.0, 0.0, 0.0};
  double denominator = 0.0;

  double du_numerator[3] = {0.0, 0.0, 0.0};
  double du_denominator = 0.0;

  double dv_numerator[3] = {0.0, 0.0, 0.0};
  double dv_denominator = 0.0;

  for (int i = 0; i <= n; ++i) {
    double Nu = Nip(i, p, u, knotU);
    double dNu = NipDerivative(i, p, u, knotU);
    for (int j = 0; j <= m; ++j) {
      double Nv = Nip(j, q, v, knotV);
      double dNv = NipDerivative(j, q, v, knotV);
      double w = weights[i][j];
      double B = Nu * Nv * w;

      // Accumuler le numérateur et le dénominateur pour S(u,v)
      for (int k = 0; k < 3; ++k) {
          numerator[k] += B * ctrlPoints[i][j][k];
      }
      denominator += B;

      // Accumuler pour la dérivée par rapport à u
      double Bu = dNu * Nv * w;
      for (int k = 0; k < 3; ++k) {
          du_numerator[k] += Bu * ctrlPoints[i][j][k];
      }
      du_denominator += Bu;

      // Accumuler pour la dérivée par rapport à v
      double Bv = Nu * dNv * w;
      for (int k = 0; k < 3; ++k) {
          dv_numerator[k] += Bv * ctrlPoints[i][j][k];
      }
      dv_denominator += Bv;
    }
  }

  // Calcul du point S(u,v)
    point[0] = numerator[0] / denominator;
    point[1] = numerator[1] / denominator;
    point[2] = numerator[2] / denominator;

    // Calcul des dérivées partielles
    du[0] = (du_numerator[0] * denominator - numerator[0] * du_denominator) / (denominator * denominator);
    du[1] = (du_numerator[1] * denominator - numerator[1] * du_denominator) / (denominator * denominator);
    du[2] = (du_numerator[2] * denominator - numerator[2] * du_denominator) / (denominator * denominator);

    dv[0] = (dv_numerator[0] * denominator - numerator[0] * dv_denominator) / (denominator * denominator);
    dv[1] = (dv_numerator[1] * denominator - numerator[1] * dv_denominator) / (denominator * denominator);
    dv[2] = (dv_numerator[2] * denominator - numerator[2] * dv_denominator) / (denominator * denominator);

    // Calcul de la normale (produit vectoriel)
    normale[0] = du[1] * dv[2] - du[2] * dv[1];
    normale[1] = du[2] * dv[0] - du[0] * dv[2];
    normale[2] = du[0] * dv[1] - du[1] * dv[0];

    // Normalisation de la normale
    float norme = sqrt(normale[0]*normale[0] + normale[1]*normale[1] + normale[2]*normale[2]);
    if (norme != 0.0f) {
        normale[0] /= norme;
        normale[1] /= norme;
        normale[2] /= norme;
    }

    // Dessiner le point et les vecteurs
    glPointSize(5.0);
    glColor3f(1.0, 0.0, 0.0); // Rouge pour le point
    glBegin(GL_POINTS);
    glVertex3fv(point);
    glEnd();

    glColor3f(0.0, 1.0, 0.0); // Vert pour du
    glBegin(GL_LINES);
    glVertex3fv(point);
    glVertex3f(point[0] + du[0], point[1] + du[1], point[2] + du[2]);
    glEnd();

    glColor3f(0.0, 0.0, 1.0); // Bleu pour dv
    glBegin(GL_LINES);
    glVertex3fv(point);
    glVertex3f(point[0] + dv[0], point[1] + dv[1], point[2] + dv[2]);
    glEnd();

    glColor3f(1.0, 1.0, 0.0); // Jaune pour la normale
    glBegin(GL_LINES);
    glVertex3fv(point);
    glVertex3f(point[0] + normale[0], point[1] + normale[1], point[2] + normale[2]);
    glEnd();
}

int main(int argc, char ** argv) {
  /* initialisation de glut et creation
     de la fenetre */
  glutInit( & argc, argv);
  glutInitDisplayMode(GLUT_RGB);
  glutInitWindowPosition(200, 200);
  glutInitWindowSize(width, height);
  glutCreateWindow("ifs");

  /* Initialisation d'OpenGL */
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glColor3f(1.0, 1.0, 1.0);
  glPointSize(1.0);

  //ifs = new Ifs();
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
void affiche_repere(void) {
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
void affichage(void) {
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
  afficheVecteursTangents(u, v);
  glPopMatrix();
  glFlush();
  glutSwapBuffers();

}

//------------------------------------------------------

//------------------------------------------------------
void clavier(unsigned char touche, int x, int y) {

  switch (touche) {
  case '+': //
    t += .1;
    if (t > 1) t = 1;
    glutPostRedisplay();
    break;
  case '-': //* ajustement du t
    t -= .1;
    if (t < 0) t = 0;
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
    if (u > 1.0f) u = 1.0f;
    glutPostRedisplay();
    break;
  case 'U':
    u -= 0.05f;
    if (u > 1.0f) u = 1.0f;
    glutPostRedisplay();
    break;
  case 'v':
    v += 0.05f;
    if (v > 1.0f) v = 1.0f;
    glutPostRedisplay();
    break;
  case 'V':
    v -= 0.05f;
    if (v > 1.0f) v = 1.0f;
    glutPostRedisplay();
    break;

  case 'q':
    exit(0);
  }

}
void mouse(int button, int state, int x, int y) {
  mouseX = x;
  mouseY = y;

  if (button == GLUT_LEFT_BUTTON) {
    if (state == GLUT_DOWN) {
      mouseLeftDown = true;
    } else if (state == GLUT_UP)
      mouseLeftDown = false;
  } else if (button == GLUT_RIGHT_BUTTON) {
    if (state == GLUT_DOWN) {
      mouseRightDown = true;
    } else if (state == GLUT_UP)
      mouseRightDown = false;
  } else if (button == GLUT_MIDDLE_BUTTON) {
    if (state == GLUT_DOWN) {
      mouseMiddleDown = true;
    } else if (state == GLUT_UP)
      mouseMiddleDown = false;
  }
}

void mouseMotion(int x, int y) {
  if (mouseLeftDown) {
    cameraAngleY += (x - mouseX);
    cameraAngleX += (y - mouseY);
    mouseX = x;
    mouseY = y;
  }
  if (mouseRightDown) {
    cameraDistance += (y - mouseY) * 0.2f;
    mouseY = y;
  }

  glutPostRedisplay();
}