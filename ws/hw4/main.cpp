#include "AMPCore.h"
#include "hw/HW4.h"
#include "ManipulatorSkeleton.h" 
#include "CSpaceSkeleton.h"
#include "HelpfulClass.h"
#include <iostream>

using namespace amp;

int main(int argc, char** argv) {
    /* Include this line to have different randomized environments every time you run your code (NOTE: this has no affect on grade()) */
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