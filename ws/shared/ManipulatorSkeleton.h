#pragma once

#include "AMPCore.h"
#include "hw/HW4.h"

class MyManipulator2D : public amp::LinkManipulator2D {
    public:
        // Default constructor 
        MyManipulator2D();
        MyManipulator2D(const std::vector<double>& link_lengths);

        //Forward kinematics
        virtual Eigen::Vector2d getJointLocation(const amp::ManipulatorState& state, uint32_t joint_index) const override;

        // Inverse kinematics
        virtual amp::ManipulatorState getConfigurationFromIK(const Eigen::Vector2d& end_effector_location) const override;
};