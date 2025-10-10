#include "AMPCore.h"

#include "hw/HW2.h"
#include "hw/HW6.h"

#include "MyAStar.h"
#include "MyCSConstructors.h"
#include "ManipulatorSkeleton.h"

using namespace amp;

// Custom heuristic for Dijkstra's (always returns 0)
class DijkstraHeuristic : public amp::SearchHeuristic {
public:
    virtual double operator()(amp::Node node) const override {
        return 0.0;
    }
};

int main(int argc, char** argv) {
    /* Include this line to have different randomized environments every time you run your code (NOTE: this has no affect on grade()) */
    amp::RNG::seed(amp::RNG::randiUnbounded());

    // 1. Wave front algorithm
    // Get the problem from HW2 Exercise 2
    Problem2D point_problem = HW2::getWorkspace1();
    // Problem2D point_problem = HW2::getWorkspace2();

    // Grid size for 0.25 cells in a 10x10 workspace, 10/0.25 = 40
    std::size_t n_cells = 40;
    
    // Construct point-agent cspace instance and wavefront algorithm
    std::shared_ptr<MyPointAgentCSConstructor> point_agent_ctor = std::make_shared<MyPointAgentCSConstructor>(n_cells);
    std::shared_ptr<WaveFrontAlgorithm> wf_algo = std::make_shared<MyWaveFrontAlgorithm>();
    PointWaveFrontAlgorithm point_algo(wf_algo, point_agent_ctor);
    Path2D path = point_algo.plan(point_problem);
    // Visualizer::makeFigure(point_problem, path);
    // Visualizer::makeFigure(*point_algo.getCSpace(), path);

    // Calculate and print path length
    double path_length = path.length();
    std::cout << "\nEx: 1" << std::endl;
    std::cout << "Number of waypoints: " << path.waypoints.size() << std::endl;
    std::cout << "Total path length: " << path_length << std::endl;

    // 2. Manipulator Planning
    MyManipulator2D manipulator;
    Problem2D manip_problem = HW6::getHW4Problem2(); 

    // Workspace start and goal 
    Eigen::Vector2d x_start(-2.0, 0.0);
    Eigen::Vector2d x_goal(2.0, 0.0);
    
    manip_problem.q_init = manipulator.getConfigurationFromIK(x_start);
    manip_problem.q_goal = manipulator.getConfigurationFromIK(x_goal);
    std::size_t manip_n_cells = 100; 
    std::shared_ptr<MyManipulatorCSConstructor> manipulator_ctor = std::make_shared<MyManipulatorCSConstructor>(manip_n_cells);
    ManipulatorWaveFrontAlgorithm manip_algo(wf_algo, manipulator_ctor);
    
    ManipulatorTrajectory2Link trajectory = manip_algo.plan(manipulator, manip_problem);
    // Visualizer::makeFigure(manip_problem, manipulator, trajectory); 
    // Visualizer::makeFigure(*manip_algo.getCSpace(), trajectory);   

    std::cout << "\nEx: 2" << std::endl;
    if (trajectory.waypoints.empty()) {
        std::cout << "No path found for the manipulator." << std::endl;
    } else {
        std::cout << "Manipulator path found with " << trajectory.waypoints.size() << " waypoints." << std::endl;
    }
    
    // 3. A* Algorithm and Dijkstra's
    ShortestPathProblem problem = HW6::getEx3SPP();
    MyAStarAlgo algo;

    // Part (a)
    std::cout << "\nEx: 3a A* " << std::endl;
    LookupSearchHeuristic heuristic_astar = HW6::getEx3Heuristic();
    MyAStarAlgo::GraphSearchResult result_astar = algo.search(problem, heuristic_astar);
    
    // Part (c)
    std::cout << "\nEx: 3c Dijkstra" << std::endl;
    DijkstraHeuristic heuristic_dijkstra;
    MyAStarAlgo::GraphSearchResult result_dijkstra = algo.search(problem, heuristic_dijkstra);
    
    // Visualizer::saveFigures();
    
    amp::HW6::grade<PointWaveFrontAlgorithm, ManipulatorWaveFrontAlgorithm, MyAStarAlgo>("AllenDevaraj.AugustinPonraj@colorado.edu", argc, argv, std::make_tuple(wf_algo, point_agent_ctor), std::make_tuple(wf_algo, manipulator_ctor), std::make_tuple());
    return 0;
}