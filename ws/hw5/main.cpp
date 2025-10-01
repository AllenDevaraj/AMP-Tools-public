// This includes all of the necessary header files in the toolbox
#include "AMPCore.h"

// Include the correct homework header
#include "hw/HW5.h"
#include "hw/HW2.h"

// Include any custom headers you created in your workspace
#include "MyGDAlgorithm.h"

using namespace amp;

// Helper function to easily create a square obstacle centered at a point
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
    /* Include this line to have different randomized environments every time you run your code */
    amp::RNG::seed(amp::RNG::randiUnbounded());

    // ==================================================================
    // ## 1. Define the specific problem with two square obstacles     ##
    // ==================================================================
    Problem2D problem;
    problem.q_init = Eigen::Vector2d(0.0, 0.0);
    problem.q_goal = Eigen::Vector2d(10.0, 0.0);
    problem.x_min = -2.0;
    problem.x_max = 12.0;
    problem.y_min = -5.0;
    problem.y_max = 5.0;
    problem.obstacles.push_back(createSquare(4.0, 1.0, 1.0));
    problem.obstacles.push_back(createSquare(7.0, -1.0, 1.0));


    // ==================================================================
    // ## 2. Instantiate the algorithm and run it on the problem     ##
    // ==================================================================
    double d_star = 10.0; // Goal influence distance
    double zetta  = 1.0;  // Attractive gain
    double Q_star = 1.5;  // Obstacle influence distance
    double eta    = 5.0;  // Repulsive gain
    MyGDAlgorithm algo(d_star, zetta, Q_star, eta);

    // Run the planner on your manually created problem
    Path2D path = algo.plan(problem);
    Visualizer::makeFigure(problem, path);


    /*
    // The code for the random environment is now commented out:
    // ==================================================================
    MyGDAlgorithm algo_rand(1.0, 1.0, 1.0, 1.0);
    Path2D path_rand;
    Problem2D prob_rand;
    bool success = HW5::generateAndCheck(algo_rand, path_rand, prob_rand);
    Visualizer::makeFigure(prob_rand, path_rand);
    */

    // Visualize your potential function on the two-square problem
    Visualizer::makeFigure(MyPotentialFunction{}, problem, 30);
    Visualizer::saveFigures();
    
    // NOTE: The grader will still run on its own set of random problems for the official test.
    // The parameters you tuned above will be passed to it.
    HW5::grade<MyGDAlgorithm>("AllenDevaraj.AugustinPonraj@colorado.edu", argc, argv, d_star, zetta, Q_star, eta);
    return 0;
}