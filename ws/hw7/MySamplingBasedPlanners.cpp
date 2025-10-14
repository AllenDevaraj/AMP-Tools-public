#include "MySamplingBasedPlanners.h"
#include <random>
#include <vector>
#include <chrono>
#include <iostream>

// Implement your PRM algorithm here
amp::Path2D MyPRM::plan(const amp::Problem2D& problem) {
    std::cerr << "\n========== MyPRM::plan() STARTED ==========\n" << std::endl;
    std::cerr << "Parameters: n=" << n_samples_ << ", r=" << radius_ << ", smooth=" << smoothing_enabled_ << std::endl;
    
    // Initialize roadmap and node locations
    roadmap = std::make_shared<amp::Graph<double>>();
    node_locations.clear();

    // Reset random seed for different results each run
    srand(time(NULL));

    std::cerr << "Starting sampling phase..." << std::endl;

    // Step 1: Add init and goal to roadmap
    node_locations[0] = problem.q_init;
    node_locations[1] = problem.q_goal;
    amp::Node node_counter = 2;

    // Step 2: Sampling Phase - generate random collision-free samples
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
    
    std::cerr << "Generated " << node_locations.size() << " nodes (including init and goal)" << std::endl;

    // Step 3: Connection Phase - connect nearby nodes
    int edges_added = 0;
    for (const auto& pair1 : node_locations) {
        for (const auto& pair2 : node_locations) {
            if (pair1.first >= pair2.first) continue; // Avoid duplicate connections and self-loops

            Eigen::Vector2d q_a = pair1.second;
            Eigen::Vector2d q_b = pair2.second;
            double distance = (q_b - q_a).norm();

            // Check if within connection radius
            if (distance <= radius_) {
                // Check if edge is collision-free
                if (!MotionPlanningHelpers::CollisionChecker::isSegmentInCollision(q_a, q_b, problem.obstacles)) {
                    // Connect in BOTH directions for undirected graph
                    roadmap->connect(pair1.first, pair2.first, distance);
                    roadmap->connect(pair2.first, pair1.first, distance);
                    edges_added++;
                }
            }
        }
    }
    
    std::cerr << "Added " << edges_added << " edges to roadmap" << std::endl;

    // Step 4: Check if init and goal are connected to anything
    int init_connections = 0;
    int goal_connections = 0;
    for (const auto& neighbor : roadmap->children(0)) {
        init_connections++;
    }
    for (const auto& neighbor : roadmap->children(1)) {
        goal_connections++;
    }
    std::cerr << "Init node (0) has " << init_connections << " connections" << std::endl;
    std::cerr << "Goal node (1) has " << goal_connections << " connections" << std::endl;

    // Step 5: Graph Search Phase
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
        std::cerr << "Original path has " << path.waypoints.size() << " waypoints" << std::endl;
    } else {
        std::cerr << "A* search failed!" << std::endl;
        return amp::Path2D(); // Return empty path if A* fails
    }

    // =========================================================
    // Path Smoothing: Shortcut Smoothing
    // =========================================================
    if (smoothing_enabled_ && path.waypoints.size() > 2) {
        std::cerr << "Starting path smoothing..." << std::endl;
        amp::Path2D smoothed_path;
        smoothed_path.waypoints.push_back(path.waypoints.front()); // Add start point

        size_t current_idx = 0;
        size_t goal_idx = path.waypoints.size() - 1;

        while (current_idx < goal_idx) {
            // Try to connect directly to the goal first
            if (!MotionPlanningHelpers::CollisionChecker::isSegmentInCollision(
                    path.waypoints[current_idx], path.waypoints[goal_idx], problem.obstacles)) {
                smoothed_path.waypoints.push_back(path.waypoints[goal_idx]);
                std::cerr << "Connected directly to goal from waypoint " << current_idx << std::endl;
                break; // We're done!
            }

            // Otherwise, find the furthest waypoint we can connect to
            size_t best_next = current_idx + 1;
            for (size_t test_idx = goal_idx - 1; test_idx > current_idx + 1; --test_idx) {
                if (!MotionPlanningHelpers::CollisionChecker::isSegmentInCollision(
                        path.waypoints[current_idx], path.waypoints[test_idx], problem.obstacles)) {
                    best_next = test_idx;
                    break;
                }
            }

            std::cerr << "Jumping from waypoint " << current_idx << " to " << best_next << std::endl;
            smoothed_path.waypoints.push_back(path.waypoints[best_next]);
            current_idx = best_next;
        }

        std::cerr << "Smoothed path has " << smoothed_path.waypoints.size() << " waypoints" << std::endl;
        return smoothed_path;
    }

    std::cerr << "Returning original path (smoothing disabled or path too short)" << std::endl;
    return path;
}

// Implement your RRT algorithm here
amp::Path2D MyRRT::plan(const amp::Problem2D& problem) {
    amp::Path2D path;
    path.waypoints.push_back(problem.q_init);
    path.waypoints.push_back(problem.q_goal);
    return path;
}