#pragma once
#pragma once

#include "Array.h"
#include "Eddy.h"
#include "region.h"
#include "ISEM1_eddy.h"

#include <climits>
#include <malloc.h>
#include <initializer_list>
#include <cassert>
#include <cmath>
#include <exception>
#include <random>
#include <iostream>
#include <map>
#include <fstream>
#include <string>

class ISEM1_region : public region {
public:
	Array<ISEM1_eddy> eddies;

	ISEM1_region(const double& _u0, const double& _dt, const double& _x_inlet, const Array<double>& _y_inlet, const Array<double>& _z_inlet, const double& rep_radius, const double& max_radius, const double& _delta) 
		: region(_u0, _dt, _x_inlet, _y_inlet, _z_inlet, max_radius, _delta) // use max radius as "base_radius" since base_radius used to determine SE region size
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

		//eps_temp.resize({ 3 }); // temporary array for epsilon values

		instantiate_eddies();
	}

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
	}

	void increment_eddy(ISEM1_eddy& _eddy) {
		// convects eddies forwards
		in_box = _eddy.convect(x_max);
		// resets eddies if eddy is past outlet of synthetic eddy region
		if (not in_box) {
			y_temp = y_rand(mt);
			z_temp = z_rand(mt);
			radius_temp = compute_radius(y_temp, z_temp); // need expression for new radius
			u_temp = velocity_fn(y_temp, z_temp);
			// evaluate new value of epsilon
			epsilon_direction();

			_eddy.reset(_eddy.position(0) - x_size, y_temp, z_temp, radius_temp, vol_sqrt, u_temp, dt, y_inlet, z_inlet, eps_temp);
		}

		double temp_x = x_inlet(0, 0); // assume inlet is planar constant x

		// determines velocity fluctuation component contribution from given eddy
		for (size_t i{ 0 }; i < _eddy.num_nodes; i++) {
			// evaluate shape function for eddy at given node
			_eddy.shape = _eddy.shape_scaling_factor;

			//_eddy.shape *= _eddy.shape_fn(x_inlet(_eddy.nodes(i, 0), _eddy.nodes(i, 1)), _eddy.position(0));
			//_eddy.shape *= _eddy.shape_fn(y_inlet(_eddy.nodes(i, 0), _eddy.nodes(i, 1)), _eddy.position(1));
			//_eddy.shape *= _eddy.shape_fn(z_inlet(_eddy.nodes(i, 0), _eddy.nodes(i, 1)), _eddy.position(2));
			_eddy.shape *= _eddy.shape_fn(temp_x, _eddy.position(0));
			_eddy.shape *= _eddy.shape_fn(_eddy.nodes_pos(i, 0), _eddy.position(1));
			_eddy.shape *= _eddy.shape_fn(_eddy.nodes_pos(i, 1), _eddy.position(2));

			u_prime(_eddy.nodes(i, 0), _eddy.nodes(i, 1)) += (a11(i) * _eddy.epsilon(0) * _eddy.shape);
			v_prime(_eddy.nodes(i, 0), _eddy.nodes(i, 1)) += (a21(i) * _eddy.epsilon(0) + a22(i) * _eddy.epsilon(1)) * _eddy.shape;
			w_prime(_eddy.nodes(i, 0), _eddy.nodes(i, 1)) += (a31(i) * _eddy.epsilon(0) + a32(i) * _eddy.epsilon(1) + a33(i) * _eddy.epsilon(2)) * _eddy.shape;

			//u_prime(_eddy.nodes(i, 0), _eddy.nodes(i, 1)) += _eddy.epsilon(0) * _eddy.shape;
			//v_prime(_eddy.nodes(i, 0), _eddy.nodes(i, 1)) += _eddy.epsilon(1) * _eddy.shape;
			//w_prime(_eddy.nodes(i, 0), _eddy.nodes(i, 1)) += _eddy.epsilon(2) * _eddy.shape;
		}
	}

	double compute_radius(const double& y, const double& z) { // determine radius for eddy at given position
		if (y < delta) {
			return std::max(0.41 * y, d_max);
		}
		else {
			return std::max(-0.082 * y + 0.41 * delta, d_max);
		}
		//return (y > 0.6 * delta) ? (- 0.59 * y + 0.64 * delta) : (0.31 * y + 0.1 * delta);
	}

	
private:
	double u_temp;

	void instantiate_eddies() {
		std::uniform_real_distribution<double> x_dist = std::uniform_real_distribution<double>(x_min, x_max);

		for (size_t i{ 0 }; i < eddies.size; i++) {
			eddies(i) = ISEM1_eddy(max_nodes);

			x_temp = x_dist(mt);
			y_temp = y_rand(mt);
			z_temp = z_rand(mt);
			radius_temp = compute_radius(y_temp, z_temp);
			u_temp = velocity_fn(y_temp, z_temp);

			// evaluate new value of epsilon
			epsilon_direction();


			eddies(i).reset(x_temp, y_temp, z_temp, radius_temp, vol_sqrt, u_temp, dt, y_inlet, z_inlet, eps_temp);
		}
	}

	const double& velocity_fn(const double& y, const double& z) {
		if (y / delta < 1) {
			return pow(y / delta, 1 / 7) * u0;
		}
		else {
			return u0;
		}
	}

public:

};



/*
class ISEM1_region {
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
	size_t max_nodes; // maximum number of nodes to be checked per iteration

	// arrays containing x, y, z coordinates of nodes on inlet plane
	Array<double> x_inlet;
	Array<double> y_inlet;
	Array<double> z_inlet;

	// array of eddies
	Array<eddy> eddies;

	Array<double> u_prime;
	Array<double> v_prime;
	Array<double> w_prime;

	Array<double> a11;
	Array<double> a21;
	Array<double> a22;
	Array<double> a31;
	Array<double> a32;
	Array<double> a33;

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

	ISEM1_region(const double& _u0, const double& _dt, const double& _x_inlet, const Array<double> _y_inlet, const Array<double> _z_inlet, const double& rep_radius, const double& max_radius, const double& _delta) {
		u0 = _u0;
		dt = _dt;
		delta = _delta;

		y_inlet = _y_inlet;
		z_inlet = _z_inlet;

		y_max = y_inlet.max();
		y_min = y_inlet.min();
		z_max = z_inlet.max();
		z_min = z_inlet.min();

		// evaluate max grid size
		d_max = std::abs(y_inlet(1, 0) - y_inlet(0, 0));
		double d_temp;
		for (size_t i{ 1 }; i < (y_inlet.shape[0] - 1); i++) {
			d_temp = std::abs(y_inlet(i+1, 0) - y_inlet(i, 0));
			//std::cout << d_temp << std::endl;
			if (d_temp > d_max) {
				d_max = d_temp;
			}
		}
		for (size_t i{ 0 }; i < z_inlet.shape[1] - 1; i++) {
			d_temp = std::abs(z_inlet(0, i+1) - z_inlet(0, i));
			//std::cout << i << ", " << d_temp << std::endl;
			if (d_temp > d_max) {
				d_max = d_temp;
			}
		}

		x_inlet.resize({ _y_inlet.shape[0], _y_inlet.shape[1] });
		x_inlet.set(_x_inlet);

		double max_rad_temp = std::max(max_radius, d_max);

		x_max = _x_inlet + max_rad_temp;
		x_min = _x_inlet - max_rad_temp;
		x_size = std::abs(x_max - x_min);


		a11.resize({ y_inlet.shape[0], y_inlet.shape[1] });
		a21.resize({ y_inlet.shape[0], y_inlet.shape[1] });
		a22.resize({ y_inlet.shape[0], y_inlet.shape[1] });
		a31.resize({ y_inlet.shape[0], y_inlet.shape[1] });
		a32.resize({ y_inlet.shape[0], y_inlet.shape[1] });
		a33.resize({ y_inlet.shape[0], y_inlet.shape[1] });

		u_prime.resize({ y_inlet.shape[0], y_inlet.shape[1] });
		v_prime.resize({ y_inlet.shape[0], y_inlet.shape[1] });
		w_prime.resize({ y_inlet.shape[0], y_inlet.shape[1] });

		vol = std::abs((x_max - x_min) * (y_max - y_min) * (z_max - z_min));
		vol_sqrt = pow(vol, 0.5);

		size_t N = trunc(vol / pow(rep_radius, 3)); // use representative radius to determine number of eddies
		eddies.resize({ N });

		vf_scaling_factor = 1 / pow(eddies.size, 0.5);

		y_rand = std::uniform_real_distribution<double>(y_min, y_max);
		z_rand = std::uniform_real_distribution<double>(z_min, z_max);
		eps_rand = std::uniform_int_distribution<int>(0, 1);
		eps_map = std::map<int, int>{ {0, -1}, {1, 1} };
		mt = std::mt19937(1729);
		
		// need method for evaluating max_nodes
		// max_nodes = ;
		max_nodes = 2000; // temporary

		instantiate_eddies();
	}

	void increment_eddies() {
		u_prime.set(0);
		v_prime.set(0);
		w_prime.set(0);

		for (size_t idx{ 0 }; idx < eddies.size; idx++) {
			increment_eddy(eddies(idx));
		}

		for (size_t i{ 0 }; i < u_prime.size; i++) { // moved scaling factor for velocity fluctuations outside to reduce number of long operations
			u_prime(i) *= vf_scaling_factor;
			v_prime(i) *= vf_scaling_factor;
			w_prime(i) *= vf_scaling_factor;
		}
	}

	void increment_eddy(eddy& _eddy) {
		// convects eddies forwards
		in_box = _eddy.convect(x_max);
		// resets eddies if eddy is past outlet of synthetic eddy region
		if (not in_box) {
			y_temp = y_rand(mt);
			z_temp = z_rand(mt);
			radius_temp = compute_radius(y_temp, z_temp); // need expression for new radius
			//_eddy.reset(x_min, y_temp, z_temp, radius_temp, vol_sqrt);
			u_temp = velocity_fn(y_temp, z_temp);
			_eddy.reset(_eddy.position(0) - x_size, y_temp, z_temp, radius_temp, vol_sqrt, u_temp, dt, y_inlet, z_inlet);
		}

		// evaluate new value of epsilon
		_eddy.epsilon_direction(eps_rand, mt, eps_map);

		// determines velocity fluctuation component contribution from given eddy
		for (size_t i{ 0 }; i < _eddy.num_nodes; i++) {
			// evaluate shape function for eddy at given node
			_eddy.shape = _eddy.shape_scaling_factor;

			_eddy.shape *= _eddy.shape_fn(x_inlet(_eddy.nodes(i, 0), _eddy.nodes(i, 1)), _eddy.position(0));
			_eddy.shape *= _eddy.shape_fn(y_inlet(_eddy.nodes(i, 0), _eddy.nodes(i, 1)), _eddy.position(1));
			_eddy.shape *= _eddy.shape_fn(z_inlet(_eddy.nodes(i, 0), _eddy.nodes(i, 1)), _eddy.position(2));

			//u_prime(i) += 1 / pow(eddies.size, 0.5) * (a11(i) * _eddy.epsilon(0) * _eddy.shape);
			//v_prime(i) += 1 / pow(eddies.size, 0.5) * (a21(i) * _eddy.epsilon(0) + a22(i) * _eddy.epsilon(1)) * _eddy.shape;
			//w_prime(i) += 1 / pow(eddies.size, 0.5) * (a31(i) * _eddy.epsilon(0) + a32(i) * _eddy.epsilon(1) + a33(i) * _eddy.epsilon(2)) * _eddy.shape;

			u_prime(_eddy.nodes(i, 0), _eddy.nodes(i, 1)) += (a11(i) * _eddy.epsilon(0) * _eddy.shape);
			v_prime(_eddy.nodes(i, 0), _eddy.nodes(i, 1)) += (a21(i) * _eddy.epsilon(0) + a22(i) * _eddy.epsilon(1)) * _eddy.shape;
			w_prime(_eddy.nodes(i, 0), _eddy.nodes(i, 1)) += (a31(i) * _eddy.epsilon(0) + a32(i) * _eddy.epsilon(1) + a33(i) * _eddy.epsilon(2)) * _eddy.shape;
		}
	}

	/*
	void increment_eddy(eddy& _eddy) {
		// convects eddies forwards
		in_box = _eddy.convect(x_max);
		// resets eddies if eddy is past outlet of synthetic eddy region
		if (not in_box) {
			y_temp = y_rand(mt);
			z_temp = z_rand(mt);
			radius_temp = compute_radius(y_temp, z_temp); // need expression for new radius
			//_eddy.reset(x_min, y_temp, z_temp, radius_temp, vol_sqrt);
			u_temp = velocity_fn(y_temp, z_temp);
			_eddy.reset(_eddy.position(0) - x_size, y_temp, z_temp, radius_temp, vol_sqrt, u_temp, dt);
		}

		// evaluate new value of epsilon
		_eddy.epsilon_direction(eps_rand, mt, eps_map);

		// determines velocity fluctuation component contribution from given eddy
		for (size_t i{ 0 }; i < y_inlet.size; i++) {
			// evaluate shape function for eddy at given node
			_eddy.shape = _eddy.shape_scaling_factor;

			_eddy.shape *= _eddy.shape_fn(x_inlet(i), _eddy.position(0));
			_eddy.shape *= _eddy.shape_fn(y_inlet(i), _eddy.position(1));
			_eddy.shape *= _eddy.shape_fn(z_inlet(i), _eddy.position(2));

			//u_prime(i) += 1 / pow(eddies.size, 0.5) * (a11(i) * _eddy.epsilon(0) * _eddy.shape);
			//v_prime(i) += 1 / pow(eddies.size, 0.5) * (a21(i) * _eddy.epsilon(0) + a22(i) * _eddy.epsilon(1)) * _eddy.shape;
			//w_prime(i) += 1 / pow(eddies.size, 0.5) * (a31(i) * _eddy.epsilon(0) + a32(i) * _eddy.epsilon(1) + a33(i) * _eddy.epsilon(2)) * _eddy.shape;

			u_prime(i) += (a11(i) * _eddy.epsilon(0) * _eddy.shape);
			v_prime(i) += (a21(i) * _eddy.epsilon(0) + a22(i) * _eddy.epsilon(1)) * _eddy.shape;
			w_prime(i) += (a31(i) * _eddy.epsilon(0) + a32(i) * _eddy.epsilon(1) + a33(i) * _eddy.epsilon(2)) * _eddy.shape;
		}
	}

	double compute_radius(const double& y, const double& z) { // determine radius for eddy at given position
		if (y < delta) {
			return std::max(0.41 * y, d_max);
		}
		else{
			return std::max(-0.082 * y + 0.41 * delta, d_max);
		}
		//return (y > 0.6 * delta) ? (- 0.59 * y + 0.64 * delta) : (0.31 * y + 0.1 * delta);
	}

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
	}

	void set_RST(Array<double> r11, Array<double> r21, Array<double> r22, Array<double> r31, Array<double> r32, Array<double> r33) {
		for (size_t i{ 0 }; i < r11.size; i++) {
			a11(i) = pow(r11(i), 0.5);
			a21(i) = r21(i) / a11(i);
			a22(i) = pow(r22(i) - a21(i) * a21(i), 0.5);
			a31(i) = r31(i) / a11(i);
			a32(i) = (r32(i) - a31(i) * a21(i)) / a22(i);
			a33(i) = pow(r33(i) - a31(i) * a31(i) - a32(i) * a32(i), 0.5);
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


private:
	double x_temp;
	double y_temp;
	double z_temp;
	double radius_temp;
	double u_temp;

	void instantiate_eddies() {
		std::uniform_real_distribution<double> x_dist = std::uniform_real_distribution<double>(x_min, x_max);

		for (size_t i{ 0 }; i < eddies.size; i++) {
			eddies(i) = eddy(max_nodes);

			x_temp = x_dist(mt);
			y_temp = y_rand(mt);
			z_temp = z_rand(mt);
			radius_temp = compute_radius(y_temp, z_temp);
			u_temp = velocity_fn(y_temp, z_temp);

			eddies(i).reset(x_temp, y_temp, z_temp, radius_temp, vol_sqrt, u_temp, dt, y_inlet, z_inlet);
			//for (size_t j{ 0 }; j < eddies(i).num_nodes; j++) {
			//	std::cout << eddies(i).nodes(j, 0) << ", " << eddies(i).nodes(j, 1) << std::endl;
			//}
		}
	}

	const double& velocity_fn(const double& y, const double& z) {
		if (y / delta < 1) {
			return pow(y / delta, 1 / 7) * u0;
		}
		else {
			return u0;
		}
	}

public:

};*/
