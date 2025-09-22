#pragma once

#include "Array.h"

#include <climits>
#include <malloc.h>
#include <initializer_list>
#include <cassert>
#include <cmath>
#include <exception>
#include <random>
#include <iostream>
#include <map>

// can add temporary variables as private attributes?

struct eddy {
public:
    Array<double> position;
    Array<int> epsilon; // -1 or 1
    double shape_scaling_factor; // scaling factor in shape function
    double radius;
    double shape;
    double increment; // streamwise velocity

    eddy() {
        position = Array<double>({ 3 });
        epsilon = Array<int>({ 3 });
    }


    // calculate epsilon
    void epsilon_direction(std::uniform_int_distribution<int>& eps_dist, std::mt19937& urng, std::map<int, int>& eps_map) {
        // 3 allocs per function call... (not due to map) -> likely due to number generation?
        // not due to number generation alone
        for (size_t i{ 0 }; i < epsilon.shape[0]; i++) { 
            epsilon(i) = eps_map[eps_dist(urng)];
        }
    }

    void reset(const double& new_x, const double& new_y, const double& new_z, const double& new_r, const double& volume_sqrt, const double& u, const double& dt) {
        position(0) = new_x;
        position(1) = new_y;
        position(2) = new_z;
        radius = new_r;
        shape_scaling_factor = volume_sqrt / pow(new_r, 1 / 3);
        increment = u * dt;
    }

    bool convect(const double& x_max) {
        position(0) += increment;
        return (position(0) < x_max) ? true : false;
        /*
        if (position(0) > x_max) {
            return true;
        }
        else {
            return false;
        }
        */
    }

    double shape_fn(const double& ref_pos, const double& eddy_pos) {
        theta = (ref_pos - eddy_pos) / radius;

        if (std::abs(theta) < 1) {
            return pow(1.5, 0.5) * (1 - std::abs(theta));
        }
        else {
            return 0;
        }

    }

private:
    double theta;
    double theta_magn;
    bool in_box;
public:

};