#include <iostream>
#include <GL/glut.h>
#include <share/matrix.h>
#include <share/utils.h>

using namespace std;

double L(double x, Vector _x, Vector _y) {
    double res = 0;
    double l;
    size_t n = _x.n;

    for (int i = 0; i < n; i++) {
        l = 1;
        for (int j = 0; j < n; j++) {
            if (j != i) l *= (x - _x[j]) / (_x[i] - _x[j]);
        }
        res += l * _y[i];
    }
    return res;
}


double testF(double x) {
    return 1.0 / (1 + 25 * x * x);
//return sin(sin(x)); // for example
}

double spline(double x, Vector _x, Vector _y) {
    size_t n = _x.n - 1;

    Matrix A = Matrix::getE(n + 2, n + 2);
    Vector b = Vector::get0(n + 2);

    for (size_t i = 2; i <= n; i++) {
        double h_i = _x[i] - _x[i - 1];
        double h_i_d1 = _x[i - 1] - _x[i - 2];
        A[i][i - 1] = h_i_d1;
        A[i][i] = 2 * (h_i_d1 + h_i);
        A[i][i + 1] = h_i;
        b[i] = 3 * ((_y[i] - _y[i - 1]) / h_i - (_y[i - 1] - _y[i - 2]) / h_i_d1);
    }

    Vector _SSCC = Utils::solveProgon(A, b);
    Vector AA = Vector::get0(n + 1);
    Vector BB = AA;
    Vector CC = AA;
    Vector DD = AA;

    for (size_t i = 1; i < n; i++) {
        double h_i = _x[i] - _x[i - 1];
        AA[i] = _y[i - 1];
        BB[i] = (_y[i] - _y[i - 1]) / h_i - h_i * (_SSCC[i + 1] + 2 * _SSCC[i]) / 3;
        CC[i] = _SSCC[i];
        DD[i] = (_SSCC[i + 1] - _SSCC[i]) / (3 * h_i);
    }

    AA[n] = _y[n - 1];
    BB[n] = (_y[n] - _y[n - 1]) / (_x[n] - _x[n - 1]) - 2 * (_x[n] - _x[n - 1]) * _SSCC[n] / 3;
    CC[n] = _SSCC[n];
    DD[n] = -_SSCC[n] / ((_x[n] - _x[n - 1]) * 3);

    double x_min = _x[0];
    double x_max = _x[n];
    size_t pos;
    if (x <= x_min) {
        pos = 1;
    } else if (x >= x_max) {
        pos = n;
    } else {
        auto v = _x._vec;
        pos = lower_bound(v.begin(), v.end(), x) - v.begin();
    }
    double xh = (x - _x[pos - 1]);

    return AA[pos] + BB[pos] * xh + CC[pos] * xh * xh + DD[pos] * xh * xh * xh;
}

/* executed when a regular key is pressed */
void keyboardDown(unsigned char key, int x, int y) {
}

/* executed when a regular key is released */
void keyboardUp(unsigned char key, int x, int y) {
}

/* executed when a special key is pressed */
/*void keyboardSpecialDown(int k, int x, int y) {

}*/

/* executed when a special key is released */
/*void keyboardSpecialUp(int k, int x, int y) {

}*/

/* reshaped window */
void reshape(int width, int height) {
    glViewport(0, 0, (GLsizei) width, (GLsizei) height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, (GLsizei) width, (GLsizei) height, 0, 100, -100);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

/* executed when button 'button' is put into state 'state' at screen position ('x', 'y') */
void mouseClick(int button, int state, int x, int y) {
}

/* executed when the mouse moves to position ('x', 'y') */
void mouseMotion(int x, int y) {
}

/* render the scene */
void draw() {
    //glEnable(GL_MULTISAMPLE);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(150, 480 - 150, 0);
    glScalef(100, -100, 1);
    glBegin(GL_LINES);
    {
        glColor3f(1.0, 0, 0);
        glVertex3f(0, 0, 0);
        glVertex3f(1, 0, 0);
        glColor3f(0, 1.0, 0);
        glVertex3f(0, 0, 0);
        glVertex3f(0, 1, 0);
        glColor3f(0, 0, 1.0);
        glVertex3f(0, 0, 0);
        glVertex3f(0, 0, 1);
    }
    glEnd();

    glColor3f(0.1, 0.1, 0.1);
    int N = 100;
    int M = 100;
    float dx = 0.1;
    float dy = 0.1;
    float x0 = -(float) N * dx / 2;
    float y0 = -(float) M * dy / 2;
    glBegin(GL_LINES);
    {
        for (int i = 0; i <= N; ++i) {
            glVertex3f((float) i * dx + x0, y0, 1.0f);
            glVertex3f((float) i * dx + x0, y0 + (float) M * dy, 1.0f);
        }

        for (int j = 0; j <= M; ++j) {
            glVertex3f(x0, y0 + (float) j * dy, 1.0f);
            glVertex3f((float) N * dx + x0, y0 + (float) j * dy, 1.0f);
        }
    }
    glEnd();
    /* render the scene here */

    //original
    glColor3f(0.7, 0.7, 0.7);
    glBegin(GL_LINE_STRIP);
    {
        for (int i = -150; i < 150; i++) {
            double x = i / 100.0;
            glVertex3f((GLfloat) x, (GLfloat) testF(x), 0.5f);
        }
    }
    glEnd();

    const int n = 10;
    Vector x_values = Vector::get0(n + 1);
    Vector y_values = Vector::get0(n + 1);


    for (int i = 0; i <= n; i++) {
        x_values[i] = i / 5.0 - 1;
        y_values[i] = testF(x_values[i]);
    }

    //interpolate points
    glPointSize(5);
    glColor3f(1, 0.5, 0.5);
    glBegin(GL_LINES);
    {
        for (int i = 0; i <= n; i++) {
            glVertex3f((GLfloat) x_values[i], (GLfloat) (y_values[i] + 0.05), 0.7f);
            glVertex3f((GLfloat) x_values[i], (GLfloat) (y_values[i] - 0.05), 0.7f);
        }
    }
    glEnd();

    //Lagrange
    glColor3f(0, 1, 1);
    glBegin(GL_LINE_STRIP);
    {
        for (int i = -150; i < 150; i++) {
            double x = i / 100.0;
            glVertex3f((GLfloat) x, (GLfloat) L(x, x_values, y_values), 0.5f);
        }
    }
    glEnd();

    //spline
    glColor3f(1, 1, 0);
    glBegin(GL_LINE_STRIP);
    {
        for (int i = -150; i < 150; i++) {
            double x = i / 100.0;
            glVertex3f((GLfloat) x, (GLfloat) spline(x, x_values, y_values), 0.5f);
        }
    }
    glEnd();

    glFlush();

    //glDisable(GL_MULTISAMPLE);
}

/* executed when program is idle */
void idle() {

    //std::cout << "tick" << app._tick <<std::endl;
    draw();
    glutSwapBuffers();
}

/* initialize OpenGL settings */
void initGL(int width, int height) {

    reshape(width, height);

    //glEnable(GL_CULL_FACE);
    // тест прозрачности, т.е. будет учитываться
    // четвертый параметр в glColor
    glEnable(GL_ALPHA_TEST);
    // тест глубины
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    // glColor будет устанавливать
    // свойства материала
    // вам не надо дополнительно
    // вызывать glMaterialfv
    glEnable(GL_COLOR_MATERIAL);
    // разрешаем освещение
    glEnable(GL_LIGHTING);
    // включаем нулевую лампу
    glEnable(GL_LIGHT0);
    // разрешаем смешение цветов
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // устанавливаем положение нулевой лампы
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glClearDepth(1.0f);

#ifdef SMOOTH
    glEnable(GL_DITHER);
//glEnable(GL_LINE_SMOOTH);
//glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_POLYGON_SMOOTH);
    glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_POINT_SMOOTH);
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
#endif
}

int main(int argc, char **argv) {

    glutInit(&argc, argv);

    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(640, 480);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("Perspective's GLUT Template");

    // register glut call backs
    glutKeyboardFunc(keyboardDown);
    glutKeyboardUpFunc(keyboardUp);
    //glutSpecialFunc(keyboardSpecialDown);
    //glutSpecialUpFunc(keyboardSpecialUp);
    glutMouseFunc(mouseClick);
    glutMotionFunc(mouseMotion);
    //glutPassiveMotionFunc(mouseMotion);
    glutReshapeFunc(reshape);
    glutDisplayFunc(draw);
    glutIdleFunc(idle);
    glutIgnoreKeyRepeat(true); // ignore keys held down

    initGL(640, 480);

    glutMainLoop();
    return 0;
}
