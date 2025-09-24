/*
#include "AMPCore.h"
#include "hw/HW4.h"
#include "ManipulatorSkeleton.h" 
#include "CSpaceSkeleton.h"
#include "HelpfulClass.h"
#include <iostream>

using namespace amp;

int main(int argc, char** argv) {
    // Include this line to have different randomized environments every time you run your code (NOTE: this has no affect on grade()) 
    amp::RNG::seed(amp::RNG::randiUnbounded());

    // 2. a) FK 
    std::cout << "Running (a): Forward Kinematics " << std::endl;
    MyManipulator2D manipulator_a({0.5, 1.0, 0.5});
    amp::ManipulatorState state_a(3);
    state_a << M_PI/6.0, M_PI/3.0, 7.0*M_PI/4.0;

    Eigen::Vector2d end_effector_a = manipulator_a.getJointLocation(state_a, 3);
    std::cout << "End-effector location: (" << end_effector_a.x() << ", " << end_effector_a.y() << ")" << std::endl;
    
    amp::Visualizer::makeFigure(manipulator_a, state_a);

    // 2. b) IK
    std::cout << "\nRunning (b): Inverse Kinematics " << std::endl;
    MyManipulator2D manipulator_b({1.0, 0.5, 1.0});
    Eigen::Vector2d target_b(2.0, 0.0);

    amp::ManipulatorState state_b = manipulator_b.getConfigurationFromIK(target_b);
    std::cout << "Target end-effector location: (" << target_b.x() << ", " << target_b.y() << ")" << std::endl;
    std::cout << "Calculated angles (rad): theta1=" << state_b[0] << ", theta2=" << state_b[1] << ", theta3=" << state_b[2] << std::endl;
    std::cout << "Stated assumption: The final link has an absolute orientation of 0 radians." << std::endl;

    amp::Visualizer::makeFigure(manipulator_b, state_b);
    // You can visualize your cspace 
    // Visualizer::makeFigure(*cspace);

    Visualizer::saveFigures();

    // Grade method
    // amp::HW4::grade<MyManipulator2D>(cspace_constructor, "nonhuman.biologic@myspace.edu", argc, argv);
    return 0;
}
*/

#include "AMPCore.h"
#include "hw/HW4.h"
#include "ManipulatorSkeleton.h"
#include "CSpaceSkeleton.h"
#include <iostream>

// Helper function to run a test case
void runTestCase(const std::string& name, const amp::Environment2D& env) {
    std::cout << "\n--- Running Test Case: " << name << " ---" << std::endl;

    MyManipulator2D manipulator({1.0, 1.0});
    MyManipulatorCSConstructor cspace_constructor(100);
    std::unique_ptr<amp::GridCSpace2D> cspace = cspace_constructor.construct(manipulator, env);

    // FIX 1: The makeFigure function for showing a manipulator in an environment needs a third argument: the state.
    // We'll show the manipulator in its "zero" angle configuration.
    amp::ManipulatorState zero_state(manipulator.nLinks());
    zero_state.setZero();
    amp::Visualizer::makeFigure(env, manipulator, zero_state);
    
    amp::Visualizer::makeFigure(*cspace);
}

int main(int argc, char** argv) {
    MyManipulatorCSConstructor cspace_constructor(100);
    
    // === Case (a): A single triangular obstacle ===
    {
        // FIX 2: Create a named std::vector for the vertices before creating the Obstacle2D object.
        std::vector<Eigen::Vector2d> vertices_a = {Eigen::Vector2d(0.25, 0.25), Eigen::Vector2d(0.0, 0.75), Eigen::Vector2d(-0.25, 0.25)};
        amp::Environment2D env_a;
        env_a.obstacles.push_back(amp::Obstacle2D(vertices_a));
        runTestCase("(a) Triangular Obstacle", env_a);
    }

    // === Case (b): Two large rectangular obstacles ===
    {
        std::vector<Eigen::Vector2d> vertices_b1 = {Eigen::Vector2d(-0.25, 1.1), Eigen::Vector2d(-0.25, 2), Eigen::Vector2d(0.25, 2), Eigen::Vector2d(0.25, 1.1)};
        std::vector<Eigen::Vector2d> vertices_b2 = {Eigen::Vector2d(-2, -2), Eigen::Vector2d(-2, -1.8), Eigen::Vector2d(2, -1.8), Eigen::Vector2d(2, -2)};
        amp::Environment2D env_b;
        env_b.obstacles.push_back(amp::Obstacle2D(vertices_b1));
        env_b.obstacles.push_back(amp::Obstacle2D(vertices_b2));
        runTestCase("(b) Two Rectangular Obstacles", env_b);
    }

    // === Case (c): Two obstacles from previous parts ===
    {
        std::vector<Eigen::Vector2d> vertices_c1 = {Eigen::Vector2d(-0.25, 1.1), Eigen::Vector2d(-0.25, 2), Eigen::Vector2d(0.25, 2), Eigen::Vector2d(0.25, 1.1)};
        std::vector<Eigen::Vector2d> vertices_c2 = {Eigen::Vector2d(-2, -0.5), Eigen::Vector2d(-2, -0.3), Eigen::Vector2d(2, -0.3), Eigen::Vector2d(2, -0.5)};
        std::vector<Eigen::Vector2d> vertices_c3 = {Eigen::Vector2d(-2, -2), Eigen::Vector2d(-2, -1.8), Eigen::Vector2d(2, -1.8), Eigen::Vector2d(2, -2)};
        amp::Environment2D env_c;
        env_c.obstacles.push_back(amp::Obstacle2D(vertices_c1));
        env_c.obstacles.push_back(amp::Obstacle2D(vertices_c2));
        env_c.obstacles.push_back(amp::Obstacle2D(vertices_c3));
        runTestCase("(c) Three Different Obstacles", env_c);
    }

    amp::Visualizer::saveFigures();

    // Make sure to replace this with your actual email address
    amp::HW4::grade<MyManipulator2D>(cspace_constructor, "your_email@colorado.edu", argc, argv);
    
    return 0;
}