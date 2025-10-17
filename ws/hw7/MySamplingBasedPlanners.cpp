#include "MySamplingBasedPlanners.h"
#include <random>
#include <vector>
#include <limits>
#include <algorithm> 
// ============================================================================
// PRM Implementation
amp::Path2D MyPRM::plan(const amp::Problem2D& problem) {
    roadmap = std::make_shared<amp::Graph<double>>();
    node_locations.clear();

    // Step 1 & 2: Sampling
    node_locations[0] = problem.q_init;
    node_locations[1] = problem.q_goal;
    amp::Node node_counter = 2;

    for (int i = 0; i < n_samples_; ++i) {
        double x_rand = ((double)rand() / RAND_MAX) * (problem.x_max - problem.x_min) + problem.x_min;
        double y_rand = ((double)rand() / RAND_MAX) * (problem.y_max - problem.y_min) + problem.y_min;
        Eigen::Vector2d q_rand(x_rand, y_rand);

        if (!MotionPlanningHelpers::CollisionChecker::isPointInCollision(q_rand, problem.obstacles)) {
            node_locations[node_counter++] = q_rand;
        }
    }

    // Step 3: Connection
    for (const auto& pair1 : node_locations) {
        for (const auto& pair2 : node_locations) {
            if (pair1.first >= pair2.first) continue;

            Eigen::Vector2d q_a = pair1.second;
            Eigen::Vector2d q_b = pair2.second;
            double distance = (q_b - q_a).norm();

            if (distance <= radius_) {
                if (!MotionPlanningHelpers::CollisionChecker::isSegmentInCollision(q_a, q_b, problem.obstacles)) {
                    roadmap->connect(pair1.first, pair2.first, distance);
                    roadmap->connect(pair2.first, pair1.first, distance);
                }
            }
        }
    }

    // Step 4: Graph Search
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
        return amp::Path2D();
    }

    // Path Smoothing
    if (smoothing_enabled_ && path.waypoints.size() > 2) {
        amp::Path2D smoothed_path;
        smoothed_path.waypoints.push_back(path.waypoints.front());

        size_t current_idx = 0;
        while (current_idx < path.waypoints.size() - 1) {
            size_t next_idx = path.waypoints.size() - 1;
            while (next_idx > current_idx + 1) {
                if (!MotionPlanningHelpers::CollisionChecker::isSegmentInCollision(
                        path.waypoints[current_idx], path.waypoints[next_idx], problem.obstacles)) {
                    break;
                }
                next_idx--;
            }
            smoothed_path.waypoints.push_back(path.waypoints[next_idx]);
            current_idx = next_idx;
        }
        return smoothed_path;
    }

    return path;
}

// RRT Implementation
amp::Node MyRRT::findNearestNode(const Eigen::Vector2d& q_rand) {
    amp::Node nearest_node = 0;
    double min_distance = std::numeric_limits<double>::max();
    
    for (const auto& [node, position] : node_locations) {
        double distance = (position - q_rand).norm();
        if (distance < min_distance) {
            min_distance = distance;
            nearest_node = node;
        }
    }
    
    return nearest_node;
}

Eigen::Vector2d MyRRT::steer(const Eigen::Vector2d& q_near, const Eigen::Vector2d& q_rand) {
    Eigen::Vector2d direction = q_rand - q_near;
    double distance = direction.norm();
    
    if (distance <= step_size_) {
        return q_rand;
    } else {
        return q_near + (direction / distance) * step_size_;
    }
}

amp::Path2D MyRRT::plan(const amp::Problem2D& problem) {
    tree = std::make_shared<amp::Graph<double>>();
    node_locations.clear();
    
    amp::Node root_node = 0;
    node_locations[root_node] = problem.q_init;
    amp::Node node_counter = 1;
    
    amp::Node goal_node = -1;
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> x_dist(problem.x_min, problem.x_max);
    std::uniform_real_distribution<> y_dist(problem.y_min, problem.y_max);
    std::uniform_real_distribution<> prob_dist(0.0, 1.0);
    
    for (int iteration = 0; iteration < max_iterations_; ++iteration) {
        Eigen::Vector2d q_rand;
        if (prob_dist(gen) < goal_bias_) {
            q_rand = problem.q_goal;
        } else {
            q_rand = Eigen::Vector2d(x_dist(gen), y_dist(gen));
        }
        
        amp::Node nearest_node = findNearestNode(q_rand);
        Eigen::Vector2d q_near = node_locations[nearest_node];
        
        Eigen::Vector2d q_new = steer(q_near, q_rand);
        
        if (MotionPlanningHelpers::CollisionChecker::isPointInCollision(q_new, problem.obstacles)) {
            continue;
        }
        
        if (MotionPlanningHelpers::CollisionChecker::isSegmentInCollision(q_near, q_new, problem.obstacles)) {
            continue;
        }
        
        amp::Node new_node = node_counter++;
        node_locations[new_node] = q_new;
        
        double edge_distance = (q_new - q_near).norm();
        tree->connect(nearest_node, new_node, edge_distance);
        
        double distance_to_goal = (q_new - problem.q_goal).norm();
        if (distance_to_goal <= epsilon_) {
            goal_node = new_node;
            break;
        }
    }
    
    amp::Path2D path;
    if (goal_node == -1) {
        return path;
    }
    
    std::vector<amp::Node> path_nodes;
    amp::Node current_node = goal_node;
    
    while (true) {
        path_nodes.push_back(current_node);
        
        if (current_node == root_node) {
            break;
        }
        
        bool found_parent = false;
        for (const auto& [node, position] : node_locations) {
            const auto& children = tree->children(node);

            if (std::find(children.begin(), children.end(), current_node) != children.end()) {
                current_node = node; // We found the parent
                found_parent = true;
                break;
            }
        }
        
        if (!found_parent) {
            return amp::Path2D();
        }
    }
    
    std::reverse(path_nodes.begin(), path_nodes.end());
    
    for (amp::Node node : path_nodes) {
        path.waypoints.push_back(node_locations[node]);
    }
    
    if ((path.waypoints.back() - problem.q_goal).norm() > 1e-6) {
        path.waypoints.push_back(problem.q_goal);
    }
    
    return path;
}