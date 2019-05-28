#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <GLFW/glfw3.h>

#include <random>

#include "double_pendulum.hpp"


void circle(const std::pair<double, double>& center, double radius) noexcept {
    const int segments = 20;

    double delta = 2 * M_PI / segments;
    double c     = cos(delta);
    double s     = sin(delta);

    std::pair<double, double> xy{radius, 0};

    glColor3d(0, 0, 0);
    glBegin(GL_LINE_LOOP);
    for (int _ = 0; _ < segments; ++_) {
        glVertex2d(xy.first + center.first, xy.second + center.second);
        glVertex2d(center.first, center.second);

        xy = {c * xy.first - s * xy.second, s * xy.first + c * xy.second};
    }
    glEnd();
}

void pendulum_loop(GLFWwindow* window) {
    std::mt19937 generator((std::random_device())());
    std::uniform_real_distribution<> rnd(0, 2 * M_PI);

    dp::state st{{rnd(generator), rnd(generator)}, {0, 0}};

    std::pair<int, int> dim;
    glfwGetWindowSize(window, &dim.first, &dim.second);

    double len = std::min(dim.first, dim.second) / 2.5;

    dp::system ss{{1, 1}, {len, len}};

    glLineWidth(4);
    
    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT);
        glBegin(GL_LINE_STRIP);
        glColor3d(1, 0, 0);
        glVertex2d(0, 0);

        std::pair<double, double> first_edge{
            ss.length.first * sin(st.theta.first),
            ss.length.first * cos(st.theta.first)
        };

        glVertex2d(
            first_edge.first / dim.first,
            -first_edge.second / dim.second
        );

        std::pair<double, double> second_edge{
            first_edge.first  + ss.length.second * sin(st.theta.second),
            first_edge.second + ss.length.second * cos(st.theta.second)
        };

        glVertex2d(
            second_edge.first / dim.first,
            -second_edge.second / dim.second
        );

        glEnd();

        circle({0, 0}, 0.01);

        circle(
            {first_edge.first / dim.first, -first_edge.second / dim.second},
            0.01
        );

        circle(
            {second_edge.first / dim.first, -second_edge.second / dim.second},
            0.01
        );

        st = dp::advance(st, ss, 0.2);

        glfwSwapBuffers(window);
        glfwPollEvents();
    }
}

int main() noexcept {
    std::pair<int, int> dim{500, 500};

    glfwInit();
    // prohibit further complications due to resizing
    glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
    // fsaa
    glfwWindowHint(GLFW_SAMPLES, 4);
    auto window = 
        glfwCreateWindow(dim.first, dim.second, "Double Pendulum", 
            nullptr, nullptr);
    glfwMakeContextCurrent(window);    

    glfwSwapInterval(1);
    glViewport(0, 0, dim.first, dim.second);
    glClearColor(1, 1, 1, 0);
    
    pendulum_loop(window);

    glfwDestroyWindow(window);
    glfwTerminate();
}