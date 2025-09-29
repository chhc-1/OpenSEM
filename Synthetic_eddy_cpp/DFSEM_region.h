#pragma once

#include "Array.h"
#include "Eddy.h"
#include "region.h"
#include "DFSEM_eddy.h"
#include <cstdlib>
#include <random>


class DFSEM_region : public region{
public:
	Array<DFSEM_eddy> eddies;
	double min_radius;
	std::normal_distribution<double> radius_dist;

	DFSEM_region(const double& _u0, const double& _dt, const double& _x_inlet, const Array<double>& _y_inlet, const Array<double>& _z_inlet, const double& _min_radius, const double& rep_radius, const double& max_radius, const double& _delta)
		: region(_u0, _dt, _x_inlet, _y_inlet, _z_inlet, max_radius, _delta)
	{
		double max_rad_temp = std::max(max_radius, d_max);

		x_max = _x_inlet + max_rad_temp;
		x_min = _x_inlet - max_rad_temp;
		x_size = std::abs(x_max - x_min);

		vol = std::abs((x_max - x_min) * (y_max - y_min) * (z_max - z_min));
		vol_sqrt = pow(vol, 0.5);

		size_t N = trunc(vol / pow(rep_radius, 3)); // use representative radius to determine number of eddies
		eddies.resize({ N });

		vf_scaling_factor = 1 / pow(eddies.size, 0.5);

		min_radius = _min_radius;

		radius_dist = std::normal_distribution<double>(1, 1);

		//std::cout << vf_scaling_factor << std::endl;

		//eps_temp.resize({ 3 }); // temporary values for epsilon

		rad_vect_temp.resize({ 3 });
		
		instantiate_eddies();
	}

	void increment_eddies() {
		u_prime.set(0);
		v_prime.set(0);
		w_prime.set(0);

		for (size_t idx{ 0 }; idx < eddies.size; idx++) {
			increment_eddy(eddies(idx));
		}

		for (size_t idx{ 0 }; idx < u_prime.size; idx++) {
			u_prime(idx) *= vf_scaling_factor;
			v_prime(idx) *= vf_scaling_factor;
			w_prime(idx) *= vf_scaling_factor;
		}
	}

	void increment_eddy(DFSEM_eddy& _eddy) {
		in_box =  _eddy.convect(x_max);
		
		if (not in_box) {
			y_temp = y_rand(mt);
			z_temp = z_rand(mt);
			compute_radius(y_temp, z_temp);

			// evaluate new value of epsilon
			epsilon_direction();

			for (size_t m{ 0 }; m < 3; m++) {
				compute_radius(y_temp, z_temp);
				rad_vect_temp(m) = radius_temp;
			}

			_eddy.reset(_eddy.position(0) - x_size, y_temp, z_temp, rad_vect_temp, vol, u0, dt, y_inlet, z_inlet, eps_temp);
		}

		double temp_x = x_inlet(0, 0); // assume inlet is planar constant x


		for (size_t i{ 0 }; i < _eddy.num_nodes; i++) {
			_eddy.calc_r(temp_x, _eddy.nodes_pos(i, 0), _eddy.nodes_pos(i, 1));
			//_eddy.shape = _eddy.shape_scaling_factor;
			//_eddy.shape *= _eddy.shape_fn(temp_x, _eddy.nodes_pos(i, 0), _eddy.nodes_pos(i, 1));

			//u_prime(_eddy.nodes(i, 0), _eddy.nodes(i, 1)) += (_eddy.rk(1) * _eddy.epsilon(2) - _eddy.rk(2) * _eddy.epsilon(1)) * _eddy.shape / _eddy.r_magn3;
			//v_prime(_eddy.nodes(i, 0), _eddy.nodes(i, 1)) += (_eddy.rk(2) * _eddy.epsilon(0) - _eddy.rk(0) * _eddy.epsilon(2)) * _eddy.shape / _eddy.r_magn3;
			//w_prime(_eddy.nodes(i, 0), _eddy.nodes(i, 1)) += (_eddy.rk(0) * _eddy.epsilon(1) - _eddy.rk(1) * _eddy.epsilon(0)) * _eddy.shape / _eddy.r_magn3;

			u_prime(_eddy.nodes(i, 0), _eddy.nodes(i, 1)) += (_eddy.rk(1) * _eddy.epsilon(2) - _eddy.rk(2) * _eddy.epsilon(1)) * _eddy.shape_fn(_eddy.radius(0));
			v_prime(_eddy.nodes(i, 0), _eddy.nodes(i, 1)) += (_eddy.rk(2) * _eddy.epsilon(0) - _eddy.rk(0) * _eddy.epsilon(2)) * _eddy.shape_fn(_eddy.radius(1));
			w_prime(_eddy.nodes(i, 0), _eddy.nodes(i, 1)) += (_eddy.rk(0) * _eddy.epsilon(1) - _eddy.rk(1) * _eddy.epsilon(0)) * _eddy.shape_fn(_eddy.radius(2));
		}
	}

	void compute_radius(const double& y, const double& z) {
		// determine mean radius from y, z position
		if (y < 0.6 * delta) {
			mean_temp_radius = 0.31 * y + 0.1 * delta;
		}
		else {
			mean_temp_radius = std::max(min_radius * 1.1, -0.59 * y + 0.64 * delta);
		}

		radius_std = pow((mean_temp_radius * 0.001) / 3, 0.5);
		radius_dist = std::normal_distribution<double>(mean_temp_radius, radius_std);

		//std::cout << radius_std << ", " << mean_temp_radius << std::endl;

		//std::cout << radius_dist(mt) << std::endl;
		radius_temp = radius_dist(mt);
		radius_temp = std::max(std::min(radius_temp, mean_temp_radius + 2 * radius_std), radius_temp - 2 * radius_std);
		

		//std::cout << radius_temp;
	}


private:
	void instantiate_eddies() {
		std::uniform_real_distribution<double> x_dist = std::uniform_real_distribution<double>(x_min, x_max);

		for (size_t i{ 0 }; i < eddies.size; i++){
			eddies(i) = DFSEM_eddy(max_nodes);

			x_temp = x_dist(mt);
			y_temp = y_rand(mt);
			z_temp = z_rand(mt);

			compute_radius(y_temp, z_temp);

			// evaluate new value of epsilon
			epsilon_direction();

			for (size_t m{ 0 }; m < 3; m++) {
				compute_radius(y_temp, z_temp);
				rad_vect_temp(m) = radius_temp;
			}

			eddies(i).reset(x_temp, y_temp, z_temp, rad_vect_temp, vol, u0, dt, y_inlet, z_inlet, eps_temp);

			eddies(i).K = 0.01;//1e-5; /// temporary - needs normalisation
		}
	}

	double mean_temp_radius;
	double radius_std;
	//double radius_temp;
	Array<double> rad_vect_temp;

public:



};



