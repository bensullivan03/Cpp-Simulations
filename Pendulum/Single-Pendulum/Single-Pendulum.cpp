#include <GL/glut.h>
#include <cmath>
#include <iostream>

using namespace std;

// Physical constants and pendulum properties
const double g = 9.81;
const double pi = M_PI;
const double l = 0.8;
const double dt = 0.01;

// Starting conditions taken as inputs in main
double theta;
double omega;

// Using RK4 - now with 2 thetas and 2 omegas
void runge_kutta(double& theta, double& omega, double dt){
    double k1_theta = omega;
    double k1_omega = -(g/l) * sin(theta);

    double k2_theta = omega + 0.5 * dt * k1_omega;
    double k2_omega = -(g/l) * sin(theta + 0.5 * dt * k1_theta);

    double k3_theta = omega + 0.5 * dt * k2_omega;
    double k3_omega = -(g/l) * sin(theta + 0.5 * dt * k2_theta);

    double k4_theta = omega + dt * k3_omega;
    double k4_omega = -(g/l) * sin(theta + dt * k3_theta);

    theta += (dt / 6.0) * (k1_theta + 2*k2_theta + 2*k3_theta + k4_theta);
    omega += (dt / 6.0) * (k1_omega + 2*k2_omega + 2*k3_omega + k4_omega);
}

// Drawing the pendulum
void drawing() {

    // Polar to Cartesian
    double x = l * sin(theta);
    double y = -l * cos(theta);

    glClear(GL_COLOR_BUFFER_BIT);
    glLoadIdentity();

    glBegin(GL_LINES);
    glVertex2f(0.0f, 0.0f);
    glVertex2f((float)x, (float)y);
    glEnd();

    glPointSize(10.0f);
    glBegin(GL_POINTS);
    glVertex2f((float)x, (float)y);
    glEnd();

    glutSwapBuffers();
}

// Timer function for animation
void timer(int = 0) {
    runge_kutta(theta, omega, dt);
    glutPostRedisplay();
    glutTimerFunc(20, timer, 0);
}

void initOpenGL() {
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glColor3f(0.0, 0.0, 0.0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(-1.5, 1.5, -1.5, 1.5);
}

int main(int argc, char** argv) {

    printf("Enter initial angle (in rad): ");
    scanf("%lf", &theta);

    printf("Enter initial angular velocity (in rad/s): ");
    scanf("%lf", &omega);

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(800, 800);
    glutCreateWindow("Single Pendulum");

    initOpenGL();

    glutDisplayFunc(drawing);
    glutTimerFunc(0, timer, 0);

    glutMainLoop();
    return 0;
}
