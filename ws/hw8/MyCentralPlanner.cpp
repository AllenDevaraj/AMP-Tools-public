#include "MyMultiAgentPlanners.h"
#include "HelpfulClass.h"
#include <Eigen/Dense>
#include <random>
#include <limits>
#include <vector>
#include <chrono>

// EXERCISE 1: CENTRALIZED PLANNER
namespace MotionPlanningHelpers {
    double pointSegmentDistance(const Eigen::Vector2d& p, const Eigen::Vector2d& a, const Eigen::Vector2d& b, Eigen::Vector2d& closest_point) {
        Eigen::Vector2d ab = b - a;
        double t = (p - a).dot(ab) / ab.squaredNorm();
        t = std::max(0.0, std::min(1.0, t));
        closest_point = a + t * ab;
        return (p - closest_point).norm();
    }
    bool isCircleInCollision(const Eigen::Vector2d& center, double radius, const amp::Obstacle2D& obstacle) {
        if (CollisionChecker::isPointInCollision(center, {obstacle})) return true;
        const auto vertices = obstacle.verticesCCW();
        for (size_t i = 0; i < vertices.size(); ++i) {
            Eigen::Vector2d p1 = vertices[i];
            Eigen::Vector2d p2 = vertices[(i + 1) % vertices.size()];
            Eigen::Vector2d closest_point;
            if (pointSegmentDistance(center, p1, p2, closest_point) < radius) return true;
        }
        return false;
    }
}
bool isConfigInCollision(const Eigen::VectorXd& q_config, const amp::MultiAgentProblem2D& problem) {
    const int num_agents = problem.agent_properties.size();
    std::vector<Eigen::Vector2d> agent_positions;
    for (int i = 0; i < num_agents; ++i) {
        agent_positions.push_back(Eigen::Vector2d(q_config(2 * i), q_config(2 * i + 1)));
    }
    for (int i = 0; i < num_agents; ++i) {
        double radius = problem.agent_properties[i].radius;
        for (const amp::Obstacle2D& obstacle : problem.obstacles) {
            if (MotionPlanningHelpers::isCircleInCollision(agent_positions[i], radius, obstacle)) return true;
        }
    }
    for (int i = 0; i < num_agents; ++i) {
        for (int j = i + 1; j < num_agents; ++j) {
            double dist_sq = (agent_positions[i] - agent_positions[j]).squaredNorm();
            double required_dist = problem.agent_properties[i].radius + problem.agent_properties[j].radius;
            if (dist_sq < required_dist * required_dist) return true;
        }
    }
    return false;
}
bool isSegmentInCollision(const Eigen::VectorXd& q_start, const Eigen::VectorXd& q_end, const amp::MultiAgentProblem2D& problem, int num_checks = 10) {
    for (int i = 1; i <= num_checks; ++i) {
        double t = static_cast<double>(i) / num_checks;
        Eigen::VectorXd q_interp = (1.0 - t) * q_start + t * q_end;
        if (isConfigInCollision(q_interp, problem)) return true;
    }
    return false;
}

amp::MultiAgentPath2D MyCentralPlanner::plan(const amp::MultiAgentProblem2D& problem) {
    auto start_time = std::chrono::high_resolution_clock::now();
    const int num_agents = problem.agent_properties.size();
    amp::MultiAgentPath2D multi_agent_path;
    multi_agent_path.agent_paths.resize(num_agents);
    if (num_agents == 0) { 
        auto end_time = std::chrono::high_resolution_clock::now();
        last_comp_time_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();
        last_tree_size = 0;
        return multi_agent_path;
    }
    const int c_space_dim = 2 * num_agents;
    const int n_iterations = 50000; 
    const double step_size = 0.5;
    const double p_goal = 0.05;
    const double epsilon = 0.25;
    Eigen::VectorXd q_init(c_space_dim), q_goal(c_space_dim);
    for (int i = 0; i < num_agents; ++i) {
        q_init.segment<2>(2 * i) = problem.agent_properties[i].q_init;
        q_goal.segment<2>(2 * i) = problem.agent_properties[i].q_goal;
    }
    std::map<amp::Node, Eigen::VectorXd> tree_nodes;
    std::map<amp::Node, amp::Node> parent_map;
    tree_nodes[0] = q_init;
    amp::Node node_counter = 1;
    std::random_device rd;
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<> goal_dist(0.0, 1.0);
    std::uniform_real_distribution<> x_dist(problem.x_min, problem.x_max);
    std::uniform_real_distribution<> y_dist(problem.y_min, problem.y_max);
    amp::Node goal_tree_node = -1;
    for (int k = 0; k < n_iterations; ++k) {
        Eigen::VectorXd q_rand(c_space_dim);
        if (goal_dist(gen) < p_goal) q_rand = q_goal;
        else {
            for (int i = 0; i < num_agents; ++i) q_rand.segment<2>(2 * i) = Eigen::Vector2d(x_dist(gen), y_dist(gen));
        }
        amp::Node nearest_node = 0;
        double min_dist_sq = std::numeric_limits<double>::max();
        for (const auto& [node, pos] : tree_nodes) {
            double dist_sq = (pos - q_rand).squaredNorm();
            if (dist_sq < min_dist_sq) { min_dist_sq = dist_sq; nearest_node = node; }
        }
        Eigen::VectorXd q_near = tree_nodes.at(nearest_node);
        Eigen::VectorXd dir = q_rand - q_near;
        double dist = dir.norm();
        Eigen::VectorXd q_new = q_near + (dist <= step_size ? dir : (dir / dist) * step_size);
        if (!isSegmentInCollision(q_near, q_new, problem)) {
            amp::Node new_node = node_counter++;
            tree_nodes[new_node] = q_new;
            parent_map[new_node] = nearest_node;
            if ((q_new - q_goal).norm() < epsilon) { goal_tree_node = new_node; break; }
        }
    }
    if (goal_tree_node != -1) {
        for (int i = 0; i < num_agents; ++i) multi_agent_path.agent_paths[i].waypoints.clear();
        std::vector<Eigen::VectorXd> path_configs;
        amp::Node current_node = goal_tree_node;
        while (true) {
            path_configs.push_back(tree_nodes.at(current_node));
            if (current_node == 0) break;
            current_node = parent_map.at(current_node);
        }
        std::reverse(path_configs.begin(), path_configs.end());
        for (const auto& config : path_configs) {
            for (int i = 0; i < num_agents; ++i) multi_agent_path.agent_paths[i].waypoints.push_back(config.segment<2>(2 * i));
        }
        for (int i = 0; i < num_agents; ++i) multi_agent_path.agent_paths[i].waypoints.push_back(q_goal.segment<2>(2 * i));
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    last_comp_time_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();
    last_tree_size = tree_nodes.size();
    return multi_agent_path;
}


// EXERCISE 2: DECENTRALIZED PLANNER
bool isDecentralizedCollision(const Eigen::Vector2d& q_new, int new_depth, double agent_radius,
                              const std::vector<amp::Obstacle2D>& static_obstacles,
                              const std::vector<amp::Path2D>& planned_paths,
                              const std::vector<double>& other_agent_radii) 
{
    for (const auto& obs : static_obstacles) {
        if (MotionPlanningHelpers::isCircleInCollision(q_new, agent_radius, obs)) {
            return true;
        }
    }
    for (size_t j = 0; j < planned_paths.size(); ++j) {
        const amp::Path2D& other_path = planned_paths[j];
        if (other_path.waypoints.empty()) continue;
        Eigen::Vector2d other_pos;
        if (new_depth < other_path.waypoints.size()) {
            other_pos = other_path.waypoints[new_depth];
        } else {
            other_pos = other_path.waypoints.back(); 
        }
        double dist_sq = (q_new - other_pos).squaredNorm();
        double required_dist = agent_radius + other_agent_radii[j];
        if (dist_sq < required_dist * required_dist) {
            return true;
        }
    }
    return false;
}

bool isDecentralizedSegmentCollision(const Eigen::Vector2d& q_start, int start_depth,
                                     const Eigen::Vector2d& q_end, int end_depth,
                                     double agent_radius, const amp::MultiAgentProblem2D& problem,
                                     const std::vector<amp::Path2D>& planned_paths,
                                     const std::vector<double>& other_agent_radii)
{
    if (isDecentralizedCollision(q_end, end_depth, agent_radius, problem.obstacles, planned_paths, other_agent_radii)) {
        return true;
    }
    if (MotionPlanningHelpers::CollisionChecker::isSegmentInCollision(q_start, q_end, problem.obstacles)) {
        return true;
    }
    return false;
}


amp::MultiAgentPath2D MyDecentralPlanner::plan(const amp::MultiAgentProblem2D& problem) {
    auto total_start_time = std::chrono::high_resolution_clock::now();
    
    const int num_agents = problem.agent_properties.size();
    amp::MultiAgentPath2D multi_agent_path;
    multi_agent_path.agent_paths.resize(num_agents);

    std::vector<amp::Path2D> planned_paths; 
    std::vector<double> planned_radii; 

    const int n_iterations = 150000;
    const double step_size = 0.5;
    const double p_goal = 0.05;
    const double epsilon = 0.25;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> goal_dist(0.0, 1.0);
    std::uniform_real_distribution<> x_dist(problem.x_min, problem.x_max);
    std::uniform_real_distribution<> y_dist(problem.y_min, problem.y_max);

    for (int i = 0; i < num_agents; ++i) {
        Eigen::Vector2d q_init = problem.agent_properties[i].q_init;
        Eigen::Vector2d q_goal = problem.agent_properties[i].q_goal;
        double agent_radius = problem.agent_properties[i].radius;

        std::map<amp::Node, std::pair<Eigen::Vector2d, int>> tree_nodes;
        std::map<amp::Node, amp::Node> parent_map;
        amp::Node node_counter = 1;
        
        tree_nodes[0] = {q_init, 0}; 
        amp::Node goal_tree_node = -1;

        for (int k = 0; k < n_iterations; ++k) {
            Eigen::Vector2d q_rand = (goal_dist(gen) < p_goal) ? q_goal : Eigen::Vector2d(x_dist(gen), y_dist(gen));

            amp::Node nearest_node = 0;
            double min_dist_sq = std::numeric_limits<double>::max();
            for (const auto& [node, data] : tree_nodes) {
                double dist_sq = (data.first - q_rand).squaredNorm();
                if (dist_sq < min_dist_sq) {
                    min_dist_sq = dist_sq;
                    nearest_node = node;
                }
            }

            Eigen::Vector2d q_near = tree_nodes.at(nearest_node).first;
            int near_depth = tree_nodes.at(nearest_node).second;
            Eigen::Vector2d dir = q_rand - q_near;
            double dist = dir.norm();
            Eigen::Vector2d q_new = q_near + (dist <= step_size ? dir : (dir / dist) * step_size);
            int new_depth = near_depth + 1;

            if (!isDecentralizedSegmentCollision(q_near, near_depth, q_new, new_depth, agent_radius, problem, planned_paths, planned_radii)) {
                amp::Node new_node = node_counter++;
                tree_nodes[new_node] = {q_new, new_depth};
                parent_map[new_node] = nearest_node;

                if ((q_new - q_goal).norm() < epsilon) {
                    goal_tree_node = new_node;
                    break;
                }
            }
        }

        if (goal_tree_node != -1) {
            amp::Path2D agent_path;
            amp::Node current_node = goal_tree_node;
            while (true) {
                agent_path.waypoints.push_back(tree_nodes.at(current_node).first);
                if (current_node == 0) break;
                current_node = parent_map.at(current_node);
            }
            std::reverse(agent_path.waypoints.begin(), agent_path.waypoints.end());
            agent_path.waypoints.push_back(q_goal);
            
            multi_agent_path.agent_paths[i] = agent_path;
            planned_paths.push_back(agent_path); 
            planned_radii.push_back(agent_radius);
        } else {
            auto total_end_time = std::chrono::high_resolution_clock::now();
            last_comp_time_ms = std::chrono::duration<double, std::milli>(total_end_time - total_start_time).count();
            return multi_agent_path; 
        }
    }

    auto total_end_time = std::chrono::high_resolution_clock::now();
    last_comp_time_ms = std::chrono::duration<double, std::milli>(total_end_time - total_start_time).count();
    return multi_agent_path;
}