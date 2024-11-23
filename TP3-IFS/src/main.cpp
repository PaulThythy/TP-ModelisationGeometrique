#include <stdlib.h>
#include <GL/glut.h>
#include "ifs.h"


void affichage(void);

void clavier(unsigned char touche,int x,int y);
void affiche_repere(void);

void mouse(int, int, int, int);
void mouseMotion(int, int);
//void reshape(int,int);


static int nbi=5 ;
static int nb_transfo=0 ;
static Ifs *ifs=NULL;

// variables globales pour OpenGL
bool mouseLeftDown;
bool mouseRightDown;
bool mouseMiddleDown;
float mouseX, mouseY;
float cameraAngleX;
float cameraAngleY;
float cameraDistance=0.;

//------------------------------------------------------
int main(int argc,char **argv)
{

  /* initialisation de glut et creation
     de la fenetre */
  glutInit(&argc,argv);
  glutInitDisplayMode(GLUT_RGB);
  glutInitWindowPosition(200,200);
  glutInitWindowSize(600,600);
  glutCreateWindow("ifs");

  /* Initialisation d'OpenGL */
  glClearColor(0.0,0.0,0.0,0.0);
  glColor3f(1.0,1.0,1.0);
  glPointSize(1.0);
	
	ifs = new Ifs();
  /* enregistrement des fonctions de rappel */
  glutDisplayFunc(affichage);
  glutKeyboardFunc(clavier);
  glutMouseFunc(mouse);
  glutMotionFunc(mouseMotion);

// initialisation des tranformations

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	 gluPerspective(45.0f,(GLfloat)200/(GLfloat)200,0.1f,100.0f);
	glMatrixMode(GL_MODELVIEW);
	//glLoadIdentity();	
//	glScalef(.7,.7,.7);
 gluLookAt(0.,0.,2., 0.,0.,0., 0.,1.,0.);
  //  glTranslatef(0.0,0.0,-5.0);

/* Entree dans la boucle principale glut */
  glutMainLoop();
  return 0;
}
//------------------------------------------------------
void affiche_repere(void)
{
  glBegin(GL_LINES);
  glColor3f(1.0,0.0,0.0);
  glVertex2f(0.,0.);
  glVertex2f(1.,0.);
  glEnd(); 

	 glBegin(GL_LINES);
  glColor3f(0.0,1.0,0.0);
  glVertex2f(0.,0.);
  glVertex2f(0.,1.);
  glEnd(); 
   glBegin(GL_LINES);
  glColor3f(0.0,0.0,1.0);
  glVertex3f(0.,0.,0.);
  glVertex3f(0.,0.,1.);
  glEnd(); 
}

//-----------------------------------------------------


//------------------------------------------------------
void affichage(void)
{
	glMatrixMode(GL_MODELVIEW);
  /* effacement de l'image avec la couleur de fond */
	glClear(GL_COLOR_BUFFER_BIT);
	glPushMatrix();
	glTranslatef(0,0,cameraDistance);
glRotatef(cameraAngleX,1.,0.,0.)	;
glRotatef(cameraAngleY,0.,1.,0.);
	affiche_repere();
	ifs->display(nbi);
glPopMatrix();
  /* on force l'affichage du resultat */
  glFlush();
}

//------------------------------------------------------


//------------------------------------------------------
void clavier(unsigned char touche,int x,int y)
{

  switch (touche)
    {
    case '+': //* affichage du carre plein 
        nbi++;
      glutPostRedisplay();
      break;
    case '-': //* affichage du carre plein 
        nbi--;
	if (nbi < 1 ) nbi=0;
      glutPostRedisplay();
      break;
    case 'f': //* affichage en mode fil de fer 
      glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
      glutPostRedisplay();
      break;
      case 'p': //* affichage du carre plein 
      glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
      glutPostRedisplay();
      break;
  case 's' : //* Affichage en mode sommets seuls 
      glPolygonMode(GL_FRONT_AND_BACK,GL_POINT);
      glutPostRedisplay();
      break;
    case 'q' : //*la touche 'q' permet de quitter le programme 
      exit(0);
    }
    
}
void mouse(int button, int state, int x, int y)
{
    mouseX = x;
    mouseY = y;

    if(button == GLUT_LEFT_BUTTON)
    {
        if(state == GLUT_DOWN)
        {
            mouseLeftDown = true;
        }
        else if(state == GLUT_UP)
            mouseLeftDown = false;
    }

    else if(button == GLUT_RIGHT_BUTTON)
    {
        if(state == GLUT_DOWN)
        {
            mouseRightDown = true;
        }
        else if(state == GLUT_UP)
            mouseRightDown = false;
    }

    else if(button == GLUT_MIDDLE_BUTTON)
    {
        if(state == GLUT_DOWN)
        {
            mouseMiddleDown = true;
        }
        else if(state == GLUT_UP)
            mouseMiddleDown = false;
    }
}


void mouseMotion(int x, int y)
{
    if(mouseLeftDown)
    {
        cameraAngleY += (x - mouseX);
        cameraAngleX += (y - mouseY);
        mouseX = x;
        mouseY = y;
    }
    if(mouseRightDown)
    {
        cameraDistance += (y - mouseY) * 0.2f;
        mouseY = y;
    }

    glutPostRedisplay();
}

    
    
