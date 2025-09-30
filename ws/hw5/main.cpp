// This includes all of the necessary header files in the toolbox
#include "AMPCore.h"

// Include the correct homework header
#include "hw/HW5.h"
#include "hw/HW2.h"

// Include any custom headers you created in your workspace
#include "MyGDAlgorithm.h"

using namespace amp;

int main(int argc, char** argv) {
    /* Include this line to have different randomized environments every time you run your code (NOTE: this has no affect on grade()) */
    amp::RNG::seed(amp::RNG::randiUnbounded());

    // // Use WO1 from Exercise 2
    // amp::Problem2D problem = HW2::getWorkspace1();

    // // Use WO2 from Exercise 2
    // amp::Problem2D problem = HW2::getWorkspace2();

    // Test your gradient descent algorithm on a random problem.
    double d_star = 10.0; // Goal influence distance
    double zetta  = 1.0;  // Attractive gain
    double Q_star = 1.5;  // Obstacle influence distance
    double eta    = 5.0;  // Repulsive gain
    MyGDAlgorithm algo(d_star, zetta, Q_star, eta);
    Path2D path;
    Problem2D prob;
    bool success = HW5::generateAndCheck(algo, path, prob);
    Visualizer::makeFigure(prob, path);

    // Visualize your potential function
    Visualizer::makeFigure(MyPotentialFunction{}, prob, 30);
    Visualizer::saveFigures();
    
    // Arguments following argv correspond to the constructor arguments of MyGDAlgorithm:
    HW5::grade<MyGDAlgorithm>("AllenDevaraj.AugustinPonraj@colorado.edu", argc, argv, algo);
    return 0;
}



/*
// This includes all of the necessary header files in the toolbox
#include "AMPCore.h"

// Include headers for both homework assignments
#include "hw/HW5.h"
#include "hw/HW2.h"

// Include your custom HW5 algorithm
#include "MyGDAlgorithm.h"

using namespace amp;

int main(int argc, char** argv) {
    // Seed the random number generator
    amp::RNG::seed(amp::RNG::randiUnbounded());

    // ==================================================================
    // ## Central Place to Define and Tune Algorithm Parameters        ##
    // ==================================================================
    // Constructor is: MyGDAlgorithm(d_star, zetta, Q_star, eta)
    double d_star = 10.0; // Goal influence distance
    double zetta  = 1.0;  // Attractive gain
    double Q_star = 1.5;  // Obstacle influence distance
    double eta    = 5.0;  // Repulsive gain
    MyGDAlgorithm algo(d_star, zetta, Q_star, eta);


    // ==================================================================
    // ## Section 1: Standard HW5 test on a random problem             ##
    // ==================================================================
    LOG("Running standard HW5 test on a random problem...");
    {
        Path2D random_path;
        Problem2D random_prob;
        bool success = HW5::generateAndCheck(algo, random_path, random_prob);
        LOG("Random Problem Test | Valid: " << (success ? "Yes!" : "No :("));
        
        Visualizer::makeFigure(random_prob, random_path);
        Visualizer::makeFigure(MyPotentialFunction{}, random_prob);
    }
    

    // ==================================================================
    // ## Section 2: Benchmark tests on HW2 workspaces                 ##
    // ==================================================================
    LOG("Running benchmark on HW2 Workspaces...");
    {
        // Test on Workspace 1
        amp::Problem2D problem1 = HW2::getWorkspace1();
        amp::Path2D path1 = algo.plan(problem1);
        bool success1 = HW2::check(path1, problem1);
        LOG("HW2 Workspace 1 | Valid: " << (success1 ? "Yes!" : "No :("));
        LOG("HW2 Workspace 1 | Length: " << path1.length());
        Visualizer::makeFigure(problem1, path1);
    }
    {
        // Test on Workspace 2
        amp::Problem2D problem2 = HW2::getWorkspace2();
        amp::Path2D path2 = algo.plan(problem2);
        bool success2 = HW2::check(path2, problem2);
        LOG("HW2 Workspace 2 | Valid: " << (success2 ? "Yes!" : "No :("));
        LOG("HW2 Workspace 2 | Length: " << path2.length());
        Visualizer::makeFigure(problem2, path2);
    }

    // ==================================================================
    // ## Section 3: Save figures and run official grader              ##
    // ==================================================================
    Visualizer::saveFigures("hw5_combined_figs");
    
    // The grader uses its own random problems and passes your tuned parameters.
    HW5::grade<MyGDAlgorithm>("nonhuman.biologic@myspace.edu", argc, argv, d_star, zetta, Q_star, eta);
    
    return 0;
}*/