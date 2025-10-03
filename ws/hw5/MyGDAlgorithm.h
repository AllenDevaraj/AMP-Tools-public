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

class MyPotentialFunction : public amp::PotentialFunction2D {
    public:
		// Potential function
        virtual double operator()(const Eigen::Vector2d& q) const override {
            return q[0] * q[0] + q[1] * q[1];
        }
		virtual Eigen::Vector2d getGradient(const Eigen::Vector2d& q) const override {
    		return Eigen::Vector2d(2.0 * q[0], 2.0 * q[1]);
		}	

};

