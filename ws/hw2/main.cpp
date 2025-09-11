// This includes all of the necessary header files in the toolbox
#include "AMPCore.h"

// Include the correct homework header
#include "hw/HW2.h"

// Include any custom headers you created in your workspace
#include "MyBugAlgorithm.h"

using namespace amp;

namespace {

amp::Obstacle2D makeRectCCW(double x1, double y1, double x2, double y2) {
    const double lx = std::min(x1, x2), rx = std::max(x1, x2);
    const double by = std::min(y1, y2), ty = std::max(y1, y2);
    std::vector<Eigen::Vector2d> V = { {lx,by}, {lx,ty}, {rx,ty}, {rx,by} };
    return amp::Obstacle2D(V);
}

amp::Problem2D Workspace1() {
    amp::Problem2D P;
    P.x_min = -1;  P.x_max = 14;
    P.y_min = -1;  P.y_max = 14;
    P.q_init = {0,0};
    P.q_goal = {10,10};

    P.obstacles.push_back(makeRectCCW(1,1, 2,5));
    P.obstacles.push_back(makeRectCCW(3,3,   4,12));
    P.obstacles.push_back(makeRectCCW(3,12, 12,13));
    P.obstacles.push_back(makeRectCCW(12,5, 13,13));
    P.obstacles.push_back(makeRectCCW(6,5,  12,6));
    return P;
}

amp::Problem2D Workspace2() {
    amp::Problem2D P;
    P.x_min = -7;  P.x_max = 36;
    P.y_min = -7;  P.y_max = 7;
    P.q_init = {0,0};
    P.q_goal = {35,0};

    P.obstacles.push_back(makeRectCCW(-6,-6, 25,-5));  // WO1
    P.obstacles.push_back(makeRectCCW(-6,5, 30, 6));   // WO2 (top bar)
    P.obstacles.push_back(makeRectCCW(-6,-5, -5, 5));  // WO3 (left wall)
    P.obstacles.push_back(makeRectCCW( 4,-5,  5, 1));  // WO4
    P.obstacles.push_back(makeRectCCW( 9, 0, 10, 5));  // WO5
    P.obstacles.push_back(makeRectCCW(14,-5, 15, 1));  // WO6
    P.obstacles.push_back(makeRectCCW(19, 0, 20, 5));  // WO7
    P.obstacles.push_back(makeRectCCW(24,-5, 25, 1));  // WO8
    P.obstacles.push_back(makeRectCCW(29, 0, 30, 5));  // WO9
    return P;
}

} // namespace

int main(int argc, char** argv) {
    amp::RNG::seed(amp::RNG::randiUnbounded());

    // ----- choose the active workspace -----
    //amp::Problem2D problem = Workspace1();
    amp::Problem2D problem = Workspace2();

    // ----- run Bug 1 -----
    MyBugAlgorithm bug1(MyBugAlgorithm::BugMode::Bug1);   // <-- qualify BugMode
    {
        amp::Path2D path1 = bug1.plan(problem);
        bool success1 = HW2::check(path1, problem);
        LOG("Bug 1 | valid: " << (success1 ? "Yes!" : "No :("));
        LOG("Bug 1 | length: " << path1.length());
        Visualizer::makeFigure(problem, path1);
    }

    // ----- run Bug 2 -----
    MyBugAlgorithm bug2(MyBugAlgorithm::BugMode::Bug2);   // <-- qualify BugMode
    {
        amp::Path2D path2 = bug2.plan(problem);
        bool success2 = HW2::check(path2, problem);
        LOG("Bug 2 | valid: " << (success2 ? "Yes!" : "No :("));
        LOG("Bug 2 | length: " << path2.length());
        Visualizer::makeFigure(problem, path2);
    }

    // Let's get crazy and generate a random environment and test your algorithm
    /*{
        amp::Path2D path; // Make empty path, problem, and collision points, as they will be created by generateAndCheck()
        amp::Problem2D random_prob; 
        std::vector<Eigen::Vector2d> collision_points;
        bool random_trial_success = HW2::generateAndCheck(algo, path, random_prob, collision_points);
        LOG("Found valid solution in random environment: " << (random_trial_success ? "Yes!" : "No :("));

        LOG("path length: " << path.length());

        // Visualize the path environment, and any collision points with obstacles
        Visualizer::makeFigure(random_prob, path, collision_points);
    }*/

    // Save both figures
    Visualizer::saveFigures(true, "hw2_figs");

    // ----- Grading -----
    // Use the templated grader since MyBugAlgorithm doesn't derive from amp::BugAlgorithm.
    // Grade Bug 1 (default choice):
    //HW2::grade<MyBugAlgorithm>("alau2568@colorado.edu.edu", argc, argv,
    //                          MyBugAlgorithm::BugMode::Bug1);

    // If you want to grade Bug 2 instead, comment the above and use:
    HW2::grade<MyBugAlgorithm>("alau2568@colorado.edu.edu", argc, argv,
                                MyBugAlgorithm::BugMode::Bug2);

    return 0;
}
