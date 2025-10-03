// This includes all of the necessary header files in the toolbox
#include "AMPCore.h"

// Include the correct homework header
#include "hw/HW5.h"
#include "hw/HW2.h"

// Include any custom headers you created in your workspace
#include "MyGDAlgorithm.h"

using namespace amp;

amp::Obstacle2D createSquare(double centerX, double centerY, double side_length) {
    double half_side = side_length / 2.0;
    std::vector<Eigen::Vector2d> vertices;
    vertices.push_back(Eigen::Vector2d(centerX - half_side, centerY - half_side));
    vertices.push_back(Eigen::Vector2d(centerX + half_side, centerY - half_side));
    vertices.push_back(Eigen::Vector2d(centerX + half_side, centerY + half_side));
    vertices.push_back(Eigen::Vector2d(centerX - half_side, centerY + half_side));
    return amp::Obstacle2D(vertices);
}

int main(int argc, char** argv) {
    amp::RNG::seed(amp::RNG::randiUnbounded());

    // Hw5 workspace
    Problem2D problem;
    problem.q_init = Eigen::Vector2d(0.0, 0.0);
    problem.q_goal = Eigen::Vector2d(10.0, 0.0);
    problem.x_min = -2.0;
    problem.x_max = 12.0;
    problem.y_min = -5.0;
    problem.y_max = 5.0;
    problem.obstacles.push_back(createSquare(4.0, 1.0, 1.0));
    problem.obstacles.push_back(createSquare(7.0, -1.0, 1.0));

    // Random
    double d_star = 5.0; // Goal influence distance
    double zetta  = 20.0;  // Attractive gain
    double Q_star = 1;  // Obstacle influence distance
    double eta    = 0.1; // Repulsive gain 
    MyGDAlgorithm algo(d_star, zetta, Q_star, eta);

    // Option 2: HW2 Workspace 1 
    Problem2D problem = HW2::getWorkspace1();

    // // Option 3: HW2 Workspace 2 
    // Problem2D problem = HW2::getWorkspace2();

    // Run the planner on your chosen problem
    Path2D path = algo.plan(problem);
    LOG("Path length: " << path.length()); 
    // bool success = HW5::generateAndCheck(algo, path, problem);
    Visualizer::makeFigure(problem, path);

    // Visualize your potential function
    Visualizer::makeFigure(MyCombinedPotential(problem, d_star, zetta, Q_star, eta), problem, 50);
    
    Visualizer::saveFigures();
    
    HW5::grade<MyGDAlgorithm>("AllenDevaraj.AugustinPonraj@colorado.edu", argc, argv, d_star, zetta, Q_star, eta);
    return 0;
}