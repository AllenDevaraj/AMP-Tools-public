#pragma once

#include "AMPCore.h"
#include "hw/HW7.h"
#include "HelpfulClass.h"
#include "MyAStar.h"

#include <random>

class MyPRM : public amp::PRM2D {
public:
    MyPRM(int n = 200, double r = 2.0, bool smooth = false) 
        : n_samples_(n), radius_(r), smoothing_enabled_(smooth) {}

    virtual amp::Path2D plan(const amp::Problem2D& problem) override;

    std::shared_ptr<amp::Graph<double>> roadmap;
    std::map<amp::Node, Eigen::Vector2d> node_locations;

private:
    int n_samples_;
    double radius_;
    bool smoothing_enabled_;
};

class MyRRT : public amp::GoalBiasRRT2D {
public:
    MyRRT(int max_iter = 5000, double step_size = 0.5, double p_goal = 0.05, double epsilon = 0.25)
        : max_iterations_(max_iter), step_size_(step_size), goal_bias_(p_goal), epsilon_(epsilon) {}

    virtual amp::Path2D plan(const amp::Problem2D& problem) override;

    std::shared_ptr<amp::Graph<double>> tree;
    std::map<amp::Node, Eigen::Vector2d> node_locations;

private:
    int max_iterations_;
    double step_size_;
    double goal_bias_;
    double epsilon_;
    
    amp::Node findNearestNode(const Eigen::Vector2d& q_rand);
    
    Eigen::Vector2d steer(const Eigen::Vector2d& q_near, const Eigen::Vector2d& q_rand);
};