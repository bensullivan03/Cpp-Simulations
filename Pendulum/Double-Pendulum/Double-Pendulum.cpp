#include <GL/glut.h>
#include <cmath>
#include <iostream>

using namespace std;

// Physical constants and pendulum properties
// l1 is 'first' string, l2 is 'second'
const double g = 9.81;
const double pi = M_PI;
const double l1 = 0.4;
const double l2 = 0.4;
const double dt = 0.0001;

// Starting conditions taken as inputs in main
double theta1;
double omega1;
double theta2;
double omega2;

// Mass of both pendulums is now relevant
double m1 = 1;
double m2 = 2;

// Derivative of theta is omega
double dtheta_dt(double theta, double omega) {
    return omega;
}

// Derivative of omega1 and omega2 - found by considering the Lagrangian
double domega1_dt(double t1, double w1, double t2, double w2) {
    double delta = t2 - t1;
    double denom = l1 * (2 * m1 + m2 - m2 * cos(2 * delta));
    return (
        -g * (2 * m1 + m2) * sin(t1)
        - m2 * g * sin(t1 - 2 * t2)
        + 2 * sin(delta) * m2 * (w2 * w2 * l2 + w1 * w1 * l1 * cos(delta))
    ) / denom;
}

double domega2_dt(double t1, double w1, double t2, double w2) {
    double delta = t2 - t1;
    double denom = l2 * (2 * m1 + m2 - m2 * cos(2 * delta));
    return (
        -2 * sin(delta) * (
            w1 * w1 * l1 * (m1 + m2)
            + g * (m1 + m2) * cos(t1)
            + w2 * w2 * l2 * m2 * cos(delta)
        )
    ) / denom;
}

// RK - now using 2 omega and 2 theta
void runge_kutta(double& theta1, double& omega1, double& theta2, double& omega2, double dt) {
    double t1 = theta1, w1 = omega1;
    double t2 = theta2, w2 = omega2;

    double k1_theta1 = dtheta_dt(t1, w1);
    double k1_theta2 = dtheta_dt(t2, w2);
    double k1_omega1 = domega1_dt(t1, w1, t2, w2);
    double k1_omega2 = domega2_dt(t1, w1, t2, w2);

    double k2_theta1 = dtheta_dt(t1 + 0.5 * dt * k1_theta1, w1 + 0.5 * dt * k1_omega1);
    double k2_theta2 = dtheta_dt(t2 + 0.5 * dt * k1_theta2, w2 + 0.5 * dt * k1_omega2);
    double k2_omega1 = domega1_dt(
        t1 + 0.5 * dt * k1_theta1,
        w1 + 0.5 * dt * k1_omega1,
        t2 + 0.5 * dt * k1_theta2,
        w2 + 0.5 * dt * k1_omega2
    );
    double k2_omega2 = domega2_dt(
        t1 + 0.5 * dt * k1_theta1,
        w1 + 0.5 * dt * k1_omega1,
        t2 + 0.5 * dt * k1_theta2,
        w2 + 0.5 * dt * k1_omega2
    );

    double k3_theta1 = dtheta_dt(t1 + 0.5 * dt * k2_theta1, w1 + 0.5 * dt * k2_omega1);
    double k3_theta2 = dtheta_dt(t2 + 0.5 * dt * k2_theta2, w2 + 0.5 * dt * k2_omega2);
    double k3_omega1 = domega1_dt(
        t1 + 0.5 * dt * k2_theta1,
        w1 + 0.5 * dt * k2_omega1,
        t2 + 0.5 * dt * k2_theta2,
        w2 + 0.5 * dt * k2_omega2
    );
    double k3_omega2 = domega2_dt(
        t1 + 0.5 * dt * k2_theta1,
        w1 + 0.5 * dt * k2_omega1,
        t2 + 0.5 * dt * k2_theta2,
        w2 + 0.5 * dt * k2_omega2
    );

    double k4_theta1 = dtheta_dt(t1 + dt * k3_theta1, w1 + dt * k3_omega1);
    double k4_theta2 = dtheta_dt(t2 + dt * k3_theta2, w2 + dt * k3_omega2);
    double k4_omega1 = domega1_dt(
        t1 + dt * k3_theta1,
        w1 + dt * k3_omega1,
        t2 + dt * k3_theta2,
        w2 + dt * k3_omega2
    );
    double k4_omega2 = domega2_dt(
        t1 + dt * k3_theta1,
        w1 + dt * k3_omega1,
        t2 + dt * k3_theta2,
        w2 + dt * k3_omega2
    );

    theta1 += (dt / 6.0) * (k1_theta1 + 2 * k2_theta1 + 2 * k3_theta1 + k4_theta1);
    omega1 += (dt / 6.0) * (k1_omega1 + 2 * k2_omega1 + 2 * k3_omega1 + k4_omega1);
    theta2 += (dt / 6.0) * (k1_theta2 + 2 * k2_theta2 + 2 * k3_theta2 + k4_theta2);
    omega2 += (dt / 6.0) * (k1_omega2 + 2 * k2_omega2 + 2 * k3_omega2 + k4_omega2);
}

// How to make the bob look like a circle
void draw_circle(float x, float y, float r, int num_segments = 50) {
    glBegin(GL_TRIANGLE_FAN);
    glVertex2f(x, y);
    for (int i = 0; i <= num_segments; ++i) {
        float angle = 2.0f * M_PI * i / num_segments;
        glVertex2f(x + r * cos(angle), y + r * sin(angle));
    }
    glEnd();
}

// Drawing the pendulum
void drawing() {
    
    // Position of first bob 
    double x1 = l1 * sin(theta1);
    double y1 = -l1 * cos(theta1);

    // Position of second bob
    double x2 = x1 + l2 * sin(theta2);
    double y2 = y1 - l2 * cos(theta2);

    glClear(GL_COLOR_BUFFER_BIT);
    glLoadIdentity();

    glBegin(GL_LINES);
    glVertex2f(0.0f, 0.0f);
    glVertex2f((float)x1, (float)y1);
    glVertex2f((float)x1, (float)y1);
    glVertex2f((float)x2, (float)y2);
    glEnd();

    // Bobs are now circles - yay!
    draw_circle((float)x1, (float)y1, 0.05);
    draw_circle((float)x2, (float)y2, 0.05);

    glutSwapBuffers();
}

// Time keeping
void timer(int = 0) {
    int steps = 20; 
    for (int i = 0; i < steps; i++) {
        runge_kutta(theta1, omega1, theta2, omega2, dt);
    }
    glutPostRedisplay();
    glutTimerFunc(2, timer, 0);
}

// Setup the window
void initOpenGL() {
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glColor3f(0.0, 0.0, 0.0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(-1.5, 1.5, -1.5, 1.5);
}

int main(int argc, char** argv) {

    printf("Enter initial angle of pendulum 1 (in rad): ");
    scanf("%lf", &theta1);

    printf("Enter initial angular velocity of pendulum 1 (in rad/s): ");
    scanf("%lf", &omega1);

    printf("Enter initial angle of pendulum 2 (in rad): ");
    scanf("%lf", &theta2);

    printf("Enter initial angular velocity of pendulum 2 (in rad/s): ");
    scanf("%lf", &omega2);

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(800, 800);
    glutCreateWindow("Double Pendulum");

    initOpenGL();

    glutDisplayFunc(drawing);
    glutTimerFunc(0, timer, 0);

    glutMainLoop();
    return 0;
}
