#include "MyGDAlgorithm.h"
#include <vector>
#include <limits>
#include "HelpfulClass.h"

amp::Path2D MyGDAlgorithm::plan(const amp::Problem2D& problem) {
    amp::Path2D path;
    path.waypoints.push_back(problem.q_init);

    // Algorithm parameters
    const double step_size = 0.1;
    const double goal_tolerance = 0.25;
    const int max_iterations = 50000;
    
    Eigen::Vector2d q_current = problem.q_init;
    int current_iteration = 0;

    while ((q_current - problem.q_goal).norm() > goal_tolerance && current_iteration < max_iterations) {
        // Attractive Gradient
        Eigen::Vector2d grad_attractive;
        double dist_to_goal = (q_current - problem.q_goal).norm();
        if (dist_to_goal <= d_star) {
            grad_attractive = zetta * (q_current - problem.q_goal);
        } else {
            grad_attractive = d_star * zetta * (q_current - problem.q_goal) / dist_to_goal;
        }

        // Repulsive Gradient
        Eigen::Vector2d grad_repulsive(0.0, 0.0);
        for (const amp::Obstacle2D& obstacle : problem.obstacles) {
            Eigen::Vector2d closest_point_on_obstacle;
            double min_dist_sq = std::numeric_limits<double>::max();
            auto vertices = obstacle.verticesCCW();
            if (vertices.empty()) continue;

            for (size_t i = 0; i < vertices.size(); ++i) {
                const Eigen::Vector2d& v1 = vertices[i];
                const Eigen::Vector2d& v2 = vertices[(i + 1) % vertices.size()];
                Eigen::Vector2d edge = v2 - v1;
                Eigen::Vector2d vec_to_q = q_current - v1;
                double t = vec_to_q.dot(edge) / edge.squaredNorm();
                t = std::max(0.0, std::min(1.0, t));
                Eigen::Vector2d projection = v1 + t * edge;
                double dist_sq = (q_current - projection).squaredNorm();
                if (dist_sq < min_dist_sq) {
                    min_dist_sq = dist_sq;
                    closest_point_on_obstacle = projection;
                }
            }
            double dist_to_obs = std::sqrt(min_dist_sq);
            
            if (dist_to_obs <= Q_star) {
                Eigen::Vector2d grad_dist = (q_current - closest_point_on_obstacle).normalized();
                grad_repulsive += eta * (1.0 / Q_star - 1.0 / dist_to_obs) * (1.0 / (dist_to_obs * dist_to_obs)) * grad_dist;
            }
        }

        // Take a Step
        Eigen::Vector2d total_force = -grad_attractive - grad_repulsive;
        if (total_force.norm() < 1e-6) break;
        Eigen::Vector2d q_next = q_current + step_size * total_force.normalized();

        // Safety Check
        if (MotionPlanningHelpers::CollisionChecker::isPointInCollision(q_next, problem.obstacles)) {
            break; 
        }
        q_current = q_next;
        MotionPlanningHelpers::append(path.waypoints, q_current);
        current_iteration++;
    }

    if ((q_current - problem.q_goal).norm() <= goal_tolerance) {
        path.waypoints.push_back(problem.q_goal);
    }
    
    return path;
}