// This includes all of the necessary header files in the toolbox
#include "AMPCore.h"

// Include the correct homework header
#include "hw/HW2.h"

// Include any custom headers you created in your workspace
#include "MyBugAlgorithm.h"

using namespace amp;

int main(int argc, char** argv) {
    // Different randomized environments each run (no effect on grade())
    amp::RNG::seed(amp::RNG::randiUnbounded());

    // Use WO1 from Exercise 2
    amp::Problem2D problem = HW2::getWorkspace1();

    // Use WO2 from Exercise 2
    // amp::Problem2D problem = HW2::getWorkspace2();

    // Make a random environment (uncomment to use)
    /*
    Random2DEnvironmentSpecification spec;
    spec.max_obstacle_region_radius = 5.0;
    spec.n_obstacles = 2;
    spec.path_clearance = 0.01;
    spec.d_sep = 0.01;
    amp::Problem2D problem = EnvironmentTools::generateRandom(spec); // Random environment
    */

    // Choose side to follow
    using Follow = MyBugAlgorithm::FollowSide;
    constexpr Follow FOLLOW = Follow::Right;     // flip these directions to test
    const char* side = (FOLLOW == Follow::Right) ? "RIGHT" : "LEFT";

    // Bug 1
    MyBugAlgorithm bug1(MyBugAlgorithm::BugMode::Bug1, FOLLOW);
    {
        amp::Path2D path1 = bug1.plan(problem);
        bool success1 = HW2::check(path1, problem);
        LOG("Bug 1 | valid: " << (success1 ? "Yes!" : "No :("));
        LOG("Bug 1 | length: " << path1.length()); 
        try { Visualizer::makeFigure(problem, path1); }
        catch (...) { LOG("Visualization error occurred (Bug 1), continuing..."); }
    }

    // Bug 2
    MyBugAlgorithm bug2(MyBugAlgorithm::BugMode::Bug2, FOLLOW);
    {
        amp::Path2D path2 = bug2.plan(problem);
        bool success2 = HW2::check(path2, problem);
        LOG("Bug 2 | valid: " << (success2 ? "Yes!" : "No :("));
        LOG("Bug 2 | length: " << path2.length());
        try { Visualizer::makeFigure(problem, path2); }
        catch (...) { LOG("Visualization error occurred (Bug 2), continuing..."); }
    }

    try {
    Visualizer::saveFigures(true, "hw2_figs");
    } catch (...) {
    LOG("Save figures error occurred, continuing...");
}
    // Use the templated grader since MyBugAlgorithm does NOT derive from amp::BugAlgorithm.
    // Grade Bug 1 by default; pass constructor args (mode, side).
    // HW2::grade<MyBugAlgorithm>("AllenDevaraj.AugustinPonraj@colorado.edu", argc, argv,
    //                           MyBugAlgorithm::BugMode::Bug1, FOLLOW);

    // To grade Bug 2 instead, comment the call above and use:
    HW2::grade<MyBugAlgorithm>("AllenDevaraj.AugustinPonraj@colorado.edu", argc, argv,
                                MyBugAlgorithm::BugMode::Bug2, FOLLOW);

    return 0;
}
