#pragma once

#include "AMPCore.h"
#include "hw/HW7.h"
#include "HelpfulClass.h" // Include your helper class
#include "MyAStar.h"       // A* is needed for the graph search phase

#include <random>

class MyPRM : public amp::PRM2D {
public:
    MyPRM(int n = 200, double r = 2.0, bool smooth = false) : n_samples_(n), radius_(r), smoothing_enabled_(smooth) {}

    virtual amp::Path2D plan(const amp::Problem2D& problem) override;

    // Public members to access the roadmap for visualization
    std::shared_ptr<amp::Graph<double>> roadmap;
    std::map<amp::Node, Eigen::Vector2d> node_locations;

private:
    int n_samples_;
    double radius_;
    bool smoothing_enabled_;
};

class MyRRT : public amp::GoalBiasRRT2D {
public:
    virtual amp::Path2D plan(const amp::Problem2D& problem) override;
};