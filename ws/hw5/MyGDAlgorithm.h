#pragma once

// This includes all of the necessary header files in the toolbox
#include "AMPCore.h"

// Include the correct homework header
#include "hw/HW5.h"

class MyGDAlgorithm : public amp::GDAlgorithm {
	public:
		// Parameters
		MyGDAlgorithm(double d_star, double zetta, double Q_star, double eta) :
			d_star(d_star),
			zetta(zetta),
			Q_star(Q_star),
			eta(eta) {}

		virtual amp::Path2D plan(const amp::Problem2D& problem) override;
	private:
		double d_star, zetta, Q_star, eta;
};

class MyCombinedPotential : public amp::PotentialFunction2D {
public:
    MyCombinedPotential(const amp::Problem2D& problem, double d_star, double zetta, double Q_star, double eta)
        : m_problem(problem), d_star(d_star), zetta(zetta), Q_star(Q_star), eta(eta) {}

    virtual double operator()(const Eigen::Vector2d& q) const override {
        return 0.0;
    }

    virtual Eigen::Vector2d getGradient(const Eigen::Vector2d& q) const override {
        
        // Attractive Gradient
        Eigen::Vector2d grad_attractive;
        double dist_to_goal = (q - m_problem.q_goal).norm();
        if (dist_to_goal <= d_star) {
            grad_attractive = zetta * (q - m_problem.q_goal);
        } else {
            grad_attractive = d_star * zetta * (q - m_problem.q_goal) / dist_to_goal;
        }

        // Repulsive Gradient
        Eigen::Vector2d grad_repulsive(0.0, 0.0);
        for (const amp::Obstacle2D& obstacle : m_problem.obstacles) {
            Eigen::Vector2d closest_point_on_obstacle;
            double min_dist_sq = std::numeric_limits<double>::max();
            auto vertices = obstacle.verticesCCW();
            if (vertices.empty()) continue;

            for (size_t i = 0; i < vertices.size(); ++i) {
                const Eigen::Vector2d& v1 = vertices[i];
                const Eigen::Vector2d& v2 = vertices[(i + 1) % vertices.size()];
                Eigen::Vector2d edge = v2 - v1;
                Eigen::Vector2d vec_to_q = q - v1;
                double t = vec_to_q.dot(edge) / edge.squaredNorm();
                t = std::max(0.0, std::min(1.0, t));
                Eigen::Vector2d projection = v1 + t * edge;
                double dist_sq = (q - projection).squaredNorm();
                if (dist_sq < min_dist_sq) {
                    min_dist_sq = dist_sq;
                    closest_point_on_obstacle = projection;
                }
            }
            double dist_to_obs = std::sqrt(min_dist_sq);
            
            if (dist_to_obs <= Q_star) {
                Eigen::Vector2d grad_dist = (q - closest_point_on_obstacle).normalized();
                grad_repulsive += eta * (1.0 / Q_star - 1.0 / dist_to_obs) * (1.0 / (dist_to_obs * dist_to_obs)) * grad_dist;
            }
        }
        
        return grad_attractive + grad_repulsive;
    }

private:
    const amp::Problem2D& m_problem;
    double d_star, zetta, Q_star, eta;
};


class MyPotentialFunction : public amp::PotentialFunction2D {
    public:
		// Potential function
        virtual double operator()(const Eigen::Vector2d& q) const override {
            return q[0] * q[0] + q[1] * q[1];
        }
		virtual Eigen::Vector2d getGradient(const Eigen::Vector2d& q) const override {
    		return Eigen::Vector2d(2.0 * q[0],  2.0 * q[1]);
		}	
};

