#include "AMPCore.h"
#include "hw/HW2.h"
#include "hw/HW5.h"
#include "MySamplingBasedPlanners.h"
#include <iostream>

using namespace amp;

// Helper function to create the square obstacles for HW5
// FIXED: This now creates a vector of vertices first, then passes it
// to the Obstacle2D constructor, which is the correct public interface.
Obstacle2D createSquare(double x, double y, double size) {
    std::vector<Eigen::Vector2d> vertices;
    vertices.push_back(Eigen::Vector2d(x - size / 2, y - size / 2));
    vertices.push_back(Eigen::Vector2d(x + size / 2, y - size / 2));
    vertices.push_back(Eigen::Vector2d(x + size / 2, y + size / 2));
    vertices.push_back(Eigen::Vector2d(x - size / 2, y + size / 2));
    return Obstacle2D(vertices);
}

int main(int argc, char** argv) {
    // =========================================================
    // Part (a): Solve HW5 Workspace with PRM
    // =========================================================
    Problem2D problem_hw5;
    problem_hw5.q_init = Eigen::Vector2d(0.0, 0.0);
    problem_hw5.q_goal = Eigen::Vector2d(10.0, 0.0);
    problem_hw5.x_min = -1.0;
    problem_hw5.x_max = 11.0;
    problem_hw5.y_min = -3.0;
    problem_hw5.y_max = 3.0;
    problem_hw5.obstacles.push_back(createSquare(4.0, 1.0, 1.0));
    problem_hw5.obstacles.push_back(createSquare(7.0, -1.0, 1.0));

    std::cout << "Solving HW5 Problem..." << std::endl;
    std::cout << "Solving HW5 Problem..." << std::endl;
    std::cout.flush();
    
    // (a).i: Plot roadmap and path for n=200, r=1
    try {
        MyPRM prm_hw5(200, 1.0, false); // n=200, r=1, no smoothing
        std::cout << "Created PRM object, calling plan()..." << std::endl;
        std::cout.flush();
        Path2D path_hw5 = prm_hw5.plan(problem_hw5);
        std::cout << "Plan returned, visualizing..." << std::endl;
        std::cout.flush();
        // Visualizer::makeFigure(problem_hw5, path_hw5, *prm_hw5.roadmap, prm_hw5.node_locations);
    } catch (const std::exception& e) {
        std::cout << "EXCEPTION in PRM: " << e.what() << std::endl;
    }
    // (a).i: Plot roadmap and path for n=200, r=1
    MyPRM prm_hw5(200, 1.0, false); // n=200, r=1, no smoothing
    Path2D path_hw5 = prm_hw5.plan(problem_hw5);
    Visualizer::makeFigure(problem_hw5, path_hw5, *prm_hw5.roadmap, prm_hw5.node_locations);

    // (a).iv: Re-evaluate with path smoothing
    MyPRM prm_hw5_smooth(200, 1.0, true); // n=200, r=1, with smoothing
    Path2D smoothed_path_hw5 = prm_hw5_smooth.plan(problem_hw5);
    Visualizer::makeFigure(problem_hw5, smoothed_path_hw5);

    // =========================================================
    // Part (b): Solve HW2 Workspaces with PRM
    // =========================================================
    Problem2D problem_w1 = HW2::getWorkspace1();
    
    std::cout << "\nSolving HW2 Workspace 1..." << std::endl;
    // (b).i: Plot roadmap and path for n=200, r=2
    MyPRM prm_w1(200, 2.0, false); // n=200, r=2, no smoothing
    Path2D path_w1 = prm_w1.plan(problem_w1);
    Visualizer::makeFigure(problem_w1, path_w1, *prm_w1.roadmap, prm_w1.node_locations);

    // (b).iv: Re-evaluate with path smoothing
    MyPRM prm_w1_smooth(200, 2.0, true); // n=200, r=2, with smoothing
    Path2D smoothed_path_w1 = prm_w1_smooth.plan(problem_w1);
    Visualizer::makeFigure(problem_w1, smoothed_path_w1);

    // FIXED: Replaced showFigures() and setFigureTitle() with the correct function
    // from your toolbox, which saves/shows all figures at once.
    Visualizer::saveFigures();

    // Grade method (can be commented out during development)
    HW7::grade<MyPRM, MyRRT>("your.email@colorado.edu", argc, argv);

    return 0;
}