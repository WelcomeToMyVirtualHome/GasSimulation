#include <GL/glut.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>

GLfloat xRotated, yRotated, zRotated;
GLdouble radius=1;
	
int reftime = 10;
double a = 0.19; 
double L = 2.3;
double *x, *y, *z;
double mm, epseps, RR, ff, L_sphere, aa, TT, tautau;
int n, N, So, Sd, Sout, Sxyz, phi, theta, angle=10;
std::fstream file;

bool running = true;
void display(void);
void reshape(int x, int y);
void Timer(int iUnused);
static void Init(void);

void getParameters(std::string paramFile){
	std::ifstream fileOut(paramFile.c_str());
	fileOut>>n>>mm>>epseps>>RR>>ff>>L_sphere>>aa>>TT>>tautau>>So>>Sd>>Sout>>Sxyz;
	fileOut.close();
	N=pow(n,3);
	x = new double[N];
	y = new double[N];
	z = new double[N];
}

void rotate(int key,int x,int y){
	switch(key)
	{
  		case GLUT_KEY_LEFT:	
			phi=(phi-angle)%360;      
			break;
		case GLUT_KEY_RIGHT:
			phi=(phi+angle)%360;
			break;
      	case GLUT_KEY_UP:
			theta=(theta-angle)%360;
			break;
      	case GLUT_KEY_DOWN:
			theta=(theta+angle)%360;
			break;
		default: 
			break;
    }
}

void keyboard(unsigned char key, int x, int y)
{
	switch(key)
	{
  		case 'p':
  			running = !running;
  			break;
		default: 
			break;
	}
}


void Timer(int iUnused)
{
	glutPostRedisplay();
	glutTimerFunc(reftime, Timer, 0);
}

void update()
{
	if(!file.eof())
		for(int i = 0; i < N; i++)
			file >> x[i] >> y[i] >> z[i];
}

void display(void)
{
	glMatrixMode(GL_MODELVIEW);
 	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
  	glTranslatef(0.0,0.0,-10.0);
    glEnable(GL_DITHER);
	glShadeModel(GL_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST); 
	glEnable(GL_POINT_SMOOTH);
	glColorMaterial(GL_FRONT, GL_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);
    glColor3f(0, 1, 1); 
    glScalef(1.0,1.0,1.0);
	glRotated(phi,0,1,0);
	glPushMatrix();
	glRotated(theta, 1, 0, 0);
	for(int i = 0; i < N; i++){
		glTranslatef(x[i], y[i], z[i]);
	
		glutSolidSphere(aa/2, 20, 20);
		glTranslatef(-x[i], -y[i], -z[i]);		
	}
	if(running)
		update();
	glutSwapBuffers();
    glColor3ub(255,255,255);            	
    glutWireSphere(L_sphere, 32, 32);
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
		std::cout<<"usage: <input_parameters> <input_datafile>"<<'\n';	
		return 0;
	}
	getParameters(argv[1]);
	file.open(argv[2], std::ios::in | std::ios::out);
	
	glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInit(&argc, argv); 
    glutInitWindowSize(700,700);
    glutCreateWindow("Solid Sphere");
	Init();
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(rotate);
	glutReshapeFunc(reshape);
	Timer(0);
    glutMainLoop();
    return 0;
}