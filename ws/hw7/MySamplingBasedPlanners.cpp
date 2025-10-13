#include "MySamplingBasedPlanners.h"
#include <random>
#include <vector>
#include <chrono>

// Implement your PRM algorithm here
amp::Path2D MyPRM::plan(const amp::Problem2D& problem) {
    // Initialize roadmap and node locations
    roadmap = std::make_shared<amp::Graph<double>>();
    node_locations.clear();

    // Reset random seed for different results each run
    srand(time(NULL));

    // Step 1 & 2: Sampling Phase
    node_locations[0] = problem.q_init;
    node_locations[1] = problem.q_goal;
    amp::Node node_counter = 2;

    for (int i = 0; i < n_samples_; ++i) {
        // Generate a random sample
        double x_rand = ((double)rand() / RAND_MAX) * (problem.x_max - problem.x_min) + problem.x_min;
        double y_rand = ((double)rand() / RAND_MAX) * (problem.y_max - problem.y_min) + problem.y_min;
        Eigen::Vector2d q_rand(x_rand, y_rand);

        // Check if collision-free
        if (!MotionPlanningHelpers::CollisionChecker::isPointInCollision(q_rand, problem.obstacles)) {
            node_locations[node_counter++] = q_rand;
        }
    }

    // Step 3: Connection Phase
    for (const auto& pair1 : node_locations) {
        for (const auto& pair2 : node_locations) {
            if (pair1.first >= pair2.first) continue;

            Eigen::Vector2d q_a = pair1.second;
            Eigen::Vector2d q_b = pair2.second;
            double distance = (q_b - q_a).norm();

            if (distance <= radius_) {
                if (!MotionPlanningHelpers::CollisionChecker::isSegmentInCollision(q_a, q_b, problem.obstacles)) {
                    roadmap->connect(pair1.first, pair2.first, distance);
                }
            }
        }
    }

    // Step 4: Graph Search Phase
    MyAStarAlgo astar;
    amp::ShortestPathProblem search_problem;
    search_problem.graph = roadmap;
    search_problem.init_node = 0;
    search_problem.goal_node = 1;

    amp::AStar::GraphSearchResult result = astar.search(search_problem, amp::SearchHeuristic());

    amp::Path2D path;
    if (result.success) {
        for (amp::Node node : result.node_path) {
            path.waypoints.push_back(node_locations[node]);
        }
    } else {
        return amp::Path2D(); // Return empty path if A* fails
    }

    // =========================================================
    // Optional Post-Processing: Path Smoothing (FINAL CORRECTED LOGIC)
    // =========================================================
    // =========================================================
    // Path Smoothing: Robust Shortcut Smoothing
    // =========================================================
    if (smoothing_enabled_ && path.waypoints.size() > 2) {
        amp::Path2D smoothed_path;
        smoothed_path.waypoints.push_back(path.waypoints.front()); // start point

        size_t current_idx = 0;
        while (current_idx < path.waypoints.size() - 1) {
            size_t best_next = current_idx + 1;

            // Try to find the furthest collision-free shortcut
            for (size_t test_idx = path.waypoints.size() - 1; test_idx > current_idx + 1; --test_idx) {
                if (!MotionPlanningHelpers::CollisionChecker::isSegmentInCollision(
                        path.waypoints[current_idx], path.waypoints[test_idx], problem.obstacles)) {
                    best_next = test_idx;
                    break;
                }
            }

            smoothed_path.waypoints.push_back(path.waypoints[best_next]);
            current_idx = best_next;

            // Safety break (avoid infinite loop in degenerate cases)
            if (current_idx >= path.waypoints.size() - 1) break;
        }

        return smoothed_path;
    }

    return path;
}

// Implement your RRT algorithm here
amp::Path2D MyRRT::plan(const amp::Problem2D& problem) {
    amp::Path2D path;
    path.waypoints.push_back(problem.q_init);
    path.waypoints.push_back(problem.q_goal);
    return path;
}