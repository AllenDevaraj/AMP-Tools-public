#include "AMPCore.h"

#include "hw/HW2.h"
#include "hw/HW6.h"

#include "MyAStar.h"
#include "MyCSConstructors.h"
#include "ManipulatorSkeleton.h"

using namespace amp;

int main(int argc, char** argv) {
    /* Include this line to have different randomized environments every time you run your code (NOTE: this has no affect on grade()) */
    amp::RNG::seed(amp::RNG::randiUnbounded());

    // ========== EXERCISE 1: Wave-front algorithm ==========
    
    // Get the problem from HW2 Exercise 2
    Problem2D point_problem = HW2::getWorkspace1();

    // Set grid size. For 0.25 cells in a 10x10 workspace, we need 10/0.25 = 40 cells.
    std::size_t n_cells = 40;
    
    // Construct point-agent cspace instance and wavefront algorithm
    std::shared_ptr<MyPointAgentCSConstructor> point_agent_ctor = std::make_shared<MyPointAgentCSConstructor>(n_cells);
    std::shared_ptr<WaveFrontAlgorithm> wf_algo = std::make_shared<MyWaveFrontAlgorithm>();
    PointWaveFrontAlgorithm point_algo(wf_algo, point_agent_ctor);
    Path2D path = point_algo.plan(point_problem);
    Visualizer::makeFigure(point_problem, path);
    Visualizer::makeFigure(*point_algo.getCSpace(), path);

    // Calculate and print path length
    double path_length = path.length();
    std::cout << "\n=== [EX1] Path Analysis ===" << std::endl;
    std::cout << "Number of waypoints: " << path.waypoints.size() << std::endl;
    std::cout << "Total path length: " << path_length << std::endl;
    std::cout << "===========================\n" << std::endl;

    // ========== EXERCISE 2: Manipulator Planning ==========
    
    // Initialize your 2-link manipulator from HW4 (default is {1.0, 1.0} links)
    MyManipulator2D manipulator;
    Problem2D manip_problem = HW6::getHW4Problem2(); // Workspace environment

    // Define the workspace start and goal for the end-effector as per the exercise
    Eigen::Vector2d x_start(-2.0, 0.0);
    Eigen::Vector2d x_goal(2.0, 0.0);
    
    // Use Inverse Kinematics to find the C-space configurations for these workspace locations
    manip_problem.q_init = manipulator.getConfigurationFromIK(x_start);
    manip_problem.q_goal = manipulator.getConfigurationFromIK(x_goal);
    
    // Create the C-space constructor and the manipulator planning algorithm
    std::size_t manip_n_cells = 100; // Using a finer grid for C-space
    std::shared_ptr<MyManipulatorCSConstructor> manipulator_ctor = std::make_shared<MyManipulatorCSConstructor>(manip_n_cells);
    ManipulatorWaveFrontAlgorithm manip_algo(wf_algo, manipulator_ctor);
    ManipulatorTrajectory2Link trajectory = manip_algo.plan(manipulator, manip_problem);
    Visualizer::makeFigure(manip_problem, manipulator, trajectory); // Workspace snapshots
    Visualizer::makeFigure(*manip_algo.getCSpace(), trajectory);   // C-space path

    std::cout << "\n=== [EX2] Path Analysis ===" << std::endl;
    if (trajectory.waypoints.empty()) {
        std::cout << "No path found for the manipulator." << std::endl;
    } else {
        std::cout << "Manipulator path found with " << trajectory.waypoints.size() << " waypoints." << std::endl;
    }
    std::cout << "===========================\n" << std::endl;
    
    // ========== EXERCISE 3: A* Algorithm ==========
    std::cout << "\n=== [EX3] A* Search Results ===" << std::endl;
    ShortestPathProblem problem = HW6::getEx3SPP();
    LookupSearchHeuristic heuristic = HW6::getEx3Heuristic();
    MyAStarAlgo algo;
    MyAStarAlgo::GraphSearchResult result = algo.search(problem, heuristic);
    std::cout << "===============================\n" << std::endl;
    // Save all generated figures
    Visualizer::saveFigures();
    
    amp::HW6::grade<PointWaveFrontAlgorithm, ManipulatorWaveFrontAlgorithm, MyAStarAlgo>("AllenDevaraj.AugustinPonraj@colorado.edu", argc, argv, std::make_tuple(wf_algo, point_agent_ctor), std::make_tuple(wf_algo, manipulator_ctor), std::make_tuple());
    return 0;
}