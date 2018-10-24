#include <GL/glut.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
using namespace std;

GLfloat xRotated, yRotated, zRotated;
GLdouble radius=1;
	
const int atomCount = 125;
int reftime = 10;
double a = 0.19; 
double L = 2.3;
double x[atomCount];
double y[atomCount];
double z[atomCount];
double xfile, yfile, zfile, mm, epseps, RR, ff, L_sphere, aa, TT, tautau;
int n, N, So, Sd, Sout, Sxyz;
fstream file;
void display(void);
void reshape(int x, int y);
void Timer(int iUnused);
static void Init(void);

void getParameters(string paramFile){
	ifstream fileOut(paramFile.c_str());
	fileOut>>n>>mm>>epseps>>RR>>ff>>L_sphere>>aa>>TT>>tautau>>So>>Sd>>Sout>>Sxyz;
	fileOut.close();
	N=pow(n,3);
}

void Timer(int iUnused)
{
	glutPostRedisplay();
	glutTimerFunc(reftime, Timer, 0);
}

void display(void)
{
    glMatrixMode(GL_MODELVIEW);
 	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
  	glTranslatef(0.0,0.0,-8.0);
   glEnable( GL_DITHER );
	glShadeModel( GL_SMOOTH );
	glHint( GL_POINT_SMOOTH_HINT, GL_NICEST ); 
	glEnable( GL_POINT_SMOOTH );
	glColorMaterial(GL_FRONT, GL_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);
    glColor3f(0.5, 0.3, 0.2); 
    glRotatef(30,1.0,1.0,1.0);
    glScalef(1.0,1.0,1.0);
		
	if(!file.eof())
	{
		for(int i = 0; i < atomCount; i++)
		{
			file >> xfile >> yfile >> zfile;
			x[i] = xfile;
			y[i] = yfile;
			z[i] = zfile;
			glTranslatef(x[i], y[i], z[i]);
			glutSolidSphere(a, 20, 20);
			glTranslatef(-x[i], -y[i], -z[i]);		
		}
	} 
	glutSwapBuffers();
    glColor3f(0.2, 0.6, 0.2); 
	glutWireSphere(20, 32, 32);
	glFlush();
    
}
static void Init(void)
{	
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_LIGHTING); 
	glEnable(GL_LIGHT0);
	glShadeModel(GL_SMOOTH);
}

void reshape(int x, int y)
{
    glMatrixMode(GL_PROJECTION);  
    glLoadIdentity(); 
    gluPerspective(39.0,(GLdouble)x/(GLdouble)y,0.6,21.0);
    glViewport(0,0,x,y);  
    display();
}

int main (int argc, char **argv)
{
	if (argc<2){
		cout<<"usage: <input_parameters> <input_datafile>"<<endl;	
		return 0;
	}
	getParameters(argv[1]);
	file.open(argv[2], ios::in | ios::out);
	glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInit(&argc, argv); 
    glutInitWindowSize(700,700);
    glutCreateWindow("Solid Sphere");
	Init();
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
	Timer(0);
    glutMainLoop();
    return 0;
}