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


// new version of eddy (stores nodes to check internally instead of sweeping entire grid each iteration)

// base version of eddy
// for sub classes of eddy: shape functions need to be implemented, reset functions need to be implemented
struct eddy {
public:
    Array<double> position;
    Array<int> epsilon; // -1 or 1
    double shape_scaling_factor; // scaling factor in shape function
    //double radius; radius should now be a member variable declared in each sub version of eddy
    double shape;
    double increment; // streamwise velocity

    Array<double> nodes_pos; // array containing actual node positions -> store to increase speed as memory is now localised?

    Array<size_t> nodes; // array containing nodes to check each iteration, fixed size - use num_nodes to control how many nodes are required
    // nodes -> num_nodes x 2 size - column 1 contains y index, column 2 contains z index
    size_t num_nodes; // number of nodes to check each iteration
    // can add array for node locations -> then all memory accessed is in same location?

    double theta_magn;

    eddy(const size_t& max_nodes) {
        position = Array<double>({ 3 });
        epsilon = Array<int>({ 3 });
        nodes.resize({ max_nodes, 2 });
        nodes_pos.resize({ max_nodes, 2 });
        num_nodes = 0;
    }

    bool convect(const double& x_max) {
        position(0) += increment;
        return (position(0) < x_max) ? true : false;
    }

    /*
    void reset(const double& new_x, const double& new_y, const double& new_z, const double& new_r, const double& volume_sqrt, const double& u, const double& dt, const Array<double>& y_plane, const Array<double>& z_plane) {
        
        // reset eddy position
        position(0) = new_x;
        position(1) = new_y;
        position(2) = new_z;
        radius = new_r;
        shape_scaling_factor = volume_sqrt / pow(new_r, 1 / 3);
        increment = u * dt;

        // insert nodes that need to be checked now
        // nodes.set(0);
        num_nodes = 0;
        double temp_dist{ 0 };
        for (size_t j{ 0 }; j < y_plane.shape[0]; j++) {
            for (size_t k{ 0 }; k < z_plane.shape[1]; k++) {
                temp_dist = pow(pow(y_plane(j, k) - new_y, 2) + pow(z_plane(j, k) - new_z, 2), 0.5);
                if (temp_dist < radius) {
                    nodes(num_nodes, 0) = j;
                    nodes(num_nodes, 1) = k;
                    num_nodes += 1;
                }
            }
        }


    }*/


    /*
    double shape_fn(const double& ref_pos, const double& eddy_pos) {
        theta = (ref_pos - eddy_pos) / radius;

        if (std::abs(theta) < 1) {
            return pow(1.5, 0.5) * (1 - std::abs(theta));
        }
        else {
            return 0;
        }

    }*/

protected:
    double theta;
    bool in_box;
public:

};


// old version
/*
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
*/
