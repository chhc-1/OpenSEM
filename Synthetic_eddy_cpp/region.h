#pragma once

#include "Eddy.h"
#include "Array.h"
#include <cstdlib>
#include <fstream>
#include <string>
#include <random>


// base version of region class
// increment_eddies, increment_eddy, radius, instantiate_eddies functions need to be implemented in sub classes
// Array of eddies also needs to be implemented in sub classes
class region {
public:
	double x_min;
	double x_max;
	double x_size;
	double y_min;
	double y_max;
	double z_min;
	double z_max;
	double vol;
	double vol_sqrt;
	double d_max; // maximum grid size
	double base_radius;
	size_t max_nodes; // maximum number of nodes an eddy needs to check



	Array<double> x_inlet;
	Array<double> y_inlet;
	Array<double> z_inlet;
	//Array<eddy> eddies;

	Array<double> u_prime;
	Array<double> v_prime;
	Array<double> w_prime;

	Array<double> a11;
	Array<double> a21;
	Array<double> a22;
	Array<double> a31;
	Array<double> a32;
	Array<double> a33;
	Array<double> tke; // average turbulent kinetic energy of flow

	double vf_scaling_factor; // velocity fluctuation scaling factor
	double delta; // BL thickness

	double u0;
	double dt;

	bool in_box;

	std::uniform_real_distribution<double> y_rand;
	std::uniform_real_distribution<double> z_rand;
	std::uniform_int_distribution<int> eps_rand;
	std::map<int, int> eps_map; // map random generated integers to -1 and 1 for direction
	std::mt19937 mt;

	Array<int> eps_temp;

	region(const double& _u0, const double& _dt, const double& _x_inlet, const Array<double>& _y_inlet, const Array<double>& _z_inlet, const double& _radius, const double& _delta) {
		instantiate_pos(_u0, _dt, _x_inlet, _y_inlet, _z_inlet, _radius, _delta);

		resize_arrays();

		calc_d_max();

		// may need to change for each SE region
		//size_t N = trunc(vol / pow(_radius, 3)); // use representative radius to determine number of eddies
		//eddies.resize({ N });

		//vf_scaling_factor = 1 / pow(eddies.size, 0.5);

		y_rand = std::uniform_real_distribution<double>(y_min, y_max);
		z_rand = std::uniform_real_distribution<double>(z_min, z_max);
		eps_rand = std::uniform_int_distribution<int>(0, 1);
		eps_map = std::map<int, int>{ {0, -1}, {1, 1} };
		mt = std::mt19937(1729);

		// need method to evaluate max nodes
		// max_nodes = ;
		//max_nodes = 2000;
		//max_nodes = 1000;//2000;
		max_nodes = 500;
		//instantiate_eddies();
	}
	/*
	void increment_eddies() {
		u_prime.set(0);
		v_prime.set(0);
		w_prime.set(0);

		for (size_t idx{ 0 }; idx < eddies.size; idx++) {
			increment_eddy(eddies(idx));
		}

		for (size_t i{ 0 }; i < u_prime.size; i++) { // moved scaling factor for velocity fluctuations outside to reduce number of long operations
			//std::cout << vf_scaling_factor;
			//std::cout << u_prime(i) << std::endl;

			u_prime(i) *= vf_scaling_factor;
			v_prime(i) *= vf_scaling_factor;
			w_prime(i) *= vf_scaling_factor;
		}
	}*/
	/*
	void increment_eddy(eddy& _eddy) {
		// convects eddies forwards
		in_box = _eddy.convect(x_max);
		// resets eddies if eddy is past outlet of synthetic eddy region
		if (not in_box) {
			y_temp = y_rand(mt);
			z_temp = z_rand(mt);
			radius_temp = compute_radius(y_temp, z_temp); // need expression for new radius
			_eddy.reset(x_min, y_temp, z_temp, radius_temp, vol_sqrt, u0, dt, y_inlet, z_inlet);
		}

		// evaluate new value of epsilon
		_eddy.epsilon_direction(eps_rand, mt, eps_map); // 3 heap allocations here...

		//std::cout << y_inlet.size << std::endl;
		// determines velocity fluctuation component contribution from given eddy
		for (size_t i{ 0 }; i < _eddy.num_nodes; i++) {
			// evaluate shape function for eddy at given node
			_eddy.shape = _eddy.shape_scaling_factor;
			//for (size_t j{ 0 }; j < 3; j++) {
			_eddy.shape *= _eddy.shape_fn(x_inlet(_eddy.nodes(i, 0), _eddy.nodes(i, 1)), _eddy.position(0));
			_eddy.shape *= _eddy.shape_fn(y_inlet(_eddy.nodes(i, 0), _eddy.nodes(i, 1)), _eddy.position(1));
			_eddy.shape *= _eddy.shape_fn(z_inlet(_eddy.nodes(i, 0), _eddy.nodes(i, 1)), _eddy.position(2));

			u_prime(_eddy.nodes(i, 0), _eddy.nodes(i, 1)) += (a11(i) * _eddy.epsilon(0) * _eddy.shape);
			v_prime(_eddy.nodes(i, 0), _eddy.nodes(i, 1)) += (a21(i) * _eddy.epsilon(0) + a22(i) * _eddy.epsilon(1)) * _eddy.shape;
			w_prime(_eddy.nodes(i, 0), _eddy.nodes(i, 1)) += (a31(i) * _eddy.epsilon(0) + a32(i) * _eddy.epsilon(1) + a33(i) * _eddy.epsilon(2)) * _eddy.shape;
		}

	}*/

	/*
	// expect compute_radius to be overwritten
	double compute_radius(const double& y, const double& z) { // determine radius for eddy at given position
		//return (y > 0.6 * delta) ? (- 0.59 * y + 0.64 * delta) : (0.31 * y + 0.1 * delta);
		return base_radius;
	}*/

	void set_HIT_RST(const double& TI) { // sets homogenous isotropic turbulence (HIT) for Reynold's Stress Tensor (RST)
		// TI is turbulence intensity


		double u_p = pow(u0 * TI, 2);
		double u_p_sqrt = pow(u_p, 0.5);

		a11.set(u_p_sqrt);
		a21.set(0);
		a22.set(u_p_sqrt);
		a31.set(0);
		a32.set(0);
		a33.set(u_p_sqrt);
		tke.set(0.5 * (u_p + u_p + u_p));
	}


	void set_RST(Array<double> r11, Array<double> r21, Array<double> r22, Array<double> r31, Array<double> r32, Array<double> r33) {
		for (size_t i{ 0 }; i < r11.size; i++) {
			a11(i) = pow(r11(i), 0.5);
			a21(i) = r21(i) / a11(i);
			a22(i) = pow(r22(i) - a21(i) * a21(i), 0.5);
			a31(i) = r31(i) / a11(i);
			a32(i) = (r32(i) - a31(i) * a21(i)) / a22(i);
			a33(i) = pow(r33(i) - a31(i) * a31(i) - a32(i) * a32(i), 0.5);
			tke(i) = 0.5 * (r11(i) + r22(i) + r33(i));
		}
	}


	void print_flucts(const size_t& ID, const std::string& file_name) {

		std::ofstream output;

		output.open(file_name + "_uprime_" + std::to_string(ID) + ".txt");
		for (size_t j{ 0 }; j < y_inlet.shape[0]; j++) {
			//std::cout << j  << "-------------" << std::endl;
			//output << j << ", ";
			for (size_t i{ 0 }; i < y_inlet.shape[1] - 1; i++) {
				//std::cout << j * y_inlet.shape[1] + i << std::endl;
				output << u_prime(j, i) << ", ";
			}
			output << u_prime(j, y_inlet.shape[1] - 1) << std::endl;
		}
		output.close();

		output.open(file_name + "_vprime_" + std::to_string(ID) + ".txt");
		for (size_t j{ 0 }; j < y_inlet.shape[0]; j++) {
			for (size_t i{ 0 }; i < y_inlet.shape[1] - 1; i++) {
				output << v_prime(j, i) << ", ";
			}
			output << v_prime(j, y_inlet.shape[1] - 1) << std::endl;
		}
		output.close();

		output.open(file_name + "_wprime_" + std::to_string(ID) + ".txt");
		for (size_t j{ 0 }; j < y_inlet.shape[0]; j++) {
			for (size_t i{ 0 }; i < y_inlet.shape[1] - 1; i++) {
				output << w_prime(j, i) << ", ";
			}
			output << w_prime(j, y_inlet.shape[1] - 1) << std::endl;
		}
		output.close();
	}

	void epsilon_direction() {
		for (size_t d{ 0 }; d < 3; d++) {
			eps_temp(d) = eps_map[eps_rand(mt)];
		}
	}

protected:
	double x_temp;
	double y_temp;
	double z_temp;
	double radius_temp;


	/*
	void instantiate_eddies() {
		std::uniform_real_distribution<double> x_dist = std::uniform_real_distribution<double>(x_min, x_max);

		for (size_t i{ 0 }; i < eddies.size; i++) {
			eddies(i) = eddy(max_nodes);

			x_temp = x_dist(mt);
			y_temp = y_rand(mt);
			z_temp = z_rand(mt);
			radius_temp = compute_radius(y_temp, z_temp);

			eddies(i).reset(x_temp, y_temp, z_temp, radius_temp, vol_sqrt, u0, dt, y_inlet, z_inlet);
		}
	}*/

	void instantiate_pos(const double& _u0, const double& _dt, const double& _x_inlet, const Array<double>& _y_inlet, const Array<double>& _z_inlet, const double& _radius, const double& _delta) {
		u0 = _u0;
		dt = _dt;
		delta = _delta;

		x_inlet.resize({ _y_inlet.shape[0], _y_inlet.shape[1] });
		x_inlet.set(_x_inlet);
		y_inlet = _y_inlet;
		z_inlet = _z_inlet;

		x_max = _x_inlet + _radius;
		x_min = _x_inlet - _radius;
		x_size = std::abs(x_max - x_min);
		y_max = y_inlet.max() + _radius;
		//y_min = y_inlet.min() - _radius;
		y_min = y_inlet.min(); // for now assume that minimum y is a wall -> no eddies should appear below this y value
		z_max = z_inlet.max() + _radius;
		z_min = z_inlet.min() - _radius;

		base_radius = _radius; // automatic scaling for eddy size; this can be changed

		vol = std::abs((x_max - x_min) * (y_max - y_min) * (z_max - z_min));
		vol_sqrt = pow(vol, 0.5);
	}

	void resize_arrays() {
		a11.resize({ y_inlet.shape[0], y_inlet.shape[1] });
		a21.resize({ y_inlet.shape[0], y_inlet.shape[1] });
		a22.resize({ y_inlet.shape[0], y_inlet.shape[1] });
		a31.resize({ y_inlet.shape[0], y_inlet.shape[1] });
		a32.resize({ y_inlet.shape[0], y_inlet.shape[1] });
		a33.resize({ y_inlet.shape[0], y_inlet.shape[1] });
		tke.resize({ y_inlet.shape[0], y_inlet.shape[1] });

		u_prime.resize({ y_inlet.shape[0], y_inlet.shape[1] });
		v_prime.resize({ y_inlet.shape[0], y_inlet.shape[1] });
		w_prime.resize({ y_inlet.shape[0], y_inlet.shape[1] });

		eps_temp.resize({ 3, 1 });
	}

	void calc_d_max() {
		d_max = std::abs(y_inlet(1, 0) - y_inlet(0, 0));
		double d_temp;
		for (size_t i{ 1 }; i < (y_inlet.shape[0] - 1); i++) {
			d_temp = std::abs(y_inlet(i + 1, 0) - y_inlet(i, 0));
			//std::cout << d_temp << std::endl;
			if (d_temp > d_max) {
				d_max = d_temp;
			}
		}
		for (size_t i{ 0 }; i < z_inlet.shape[1] - 1; i++) {
			d_temp = std::abs(z_inlet(0, i + 1) - z_inlet(0, i));
			//std::cout << i << ", " << d_temp << std::endl;
			if (d_temp > d_max) {
				d_max = d_temp;
			}
		}
	}

public:


};