#include <GL/glut.h>
#include <vector>
#include <cmath>

// Constants of the simulation
const int N = 300;
const float c = 1.0;          
const float dx = 0.01;
const float dt = 0.001;


float u_prev[N], u_curr[N], u_next[N];

// Initialising the wave
void initwave() {
    for (int i = 0; i < N; ++i) {
        float x = (float)i / (N - 1);
        u_prev[i] = u_curr[i] = std::sin(2 * M_PI * x);
    }
    u_prev[0] = u_curr[0] = 0.0;
    u_prev[N - 1] = u_curr[N - 1] = 0.0;
}

// Updating the wave position
void updateWave() {
    for (int i = 1; i < N - 1; ++i) {
        u_next[i] = 2 * u_curr[i] - u_prev[i] + c * c * dt * dt / (dx * dx) *
                    (u_curr[i + 1] - 2 * u_curr[i] + u_curr[i - 1]);
    }
    u_next[0] = u_next[N - 1] = 0.0;

    for (int i = 0; i < N; ++i) {
        u_prev[i] = u_curr[i];
        u_curr[i] = u_next[i];
    }
}

// Drawing the wave
void drawing() {
    // Clear the screen
    glClear(GL_COLOR_BUFFER_BIT);

    // Starts drawing
    glBegin(GL_LINE_STRIP);
    for (int i = 0; i < N; ++i) {
        float x = -1.0 + 2.0 * i / (N - 1);
        // Specifies a point to draw
        glVertex2f(x, u_curr[i]);
    }
    glEnd();

    // Stops tearing
    glutSwapBuffers();
}

// Handling timekeeping
void timer(int = 0) {
    int steps = 20; 
    for (int i = 0; i < steps; i++) {
        updateWave();
    }
    glutPostRedisplay();
    // Calls timer every 10ms
    glutTimerFunc(10, timer, 0);
}

void initGL() {
    // Background colour
    glClearColor(1, 1, 1, 1);
    // Colour to draw with
    glColor3f(0.0, 0.0, 0.0);
    // Tells openGL that this is a 2D render
    glMatrixMode(GL_PROJECTION);
    // Don't understand but it seems to be needed
    glLoadIdentity();
    // Sets coordinate system range
    gluOrtho2D(-1.0, 1.0, -1.0, 1.0);
}

int main(int argc, char** argv) {
    // Initialises GLUT
    glutInit(&argc, argv);
    // Sets display properties
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    // Sets window size
    glutInitWindowSize(800, 800);
    // Creates window
    glutCreateWindow("Wave on a String");

    initGL();
    initwave();

    // Calls the drawign function
    glutDisplayFunc(drawing);

    // Continuously calling the timer function
    glutTimerFunc(0, timer, 0);

    // Keep the program going
    glutMainLoop();
    return 0;
}
