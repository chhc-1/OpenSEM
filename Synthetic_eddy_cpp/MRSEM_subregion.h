#pragma once

#include "MRSEM_eddy.h"
#include "region.h"
#include "Array.h"
#include "interpolate.h"

#include <cstdlib>

class MRSEM_subregion : public region {
public:
	Array<MRSEM_eddy> eddies;
	shape_fn shape_t[3];
	shape_fn shape_y[3];
	shape_fn shape_z[3];

	double lx;
	double lt;
	double ly;
	double lz;

	double t;

	std::uniform_real_distribution<double> t_rand;

	// can add distribution to eddy radii?

	MRSEM_subregion()
		:region()
	{

	}

	MRSEM_subregion(const double& _u0, const double& _dt, const double& _x_inlet, const Array<double>& _y_inlet, const Array<double>& _z_inlet
		, const double& _radius, const double& _delta, const double& _lx, const double& _ly, const double& _lz, const double& _y_min,
		const double& _y_max)
		: region(_u0, _dt, _x_inlet, _y_inlet, _z_inlet, _radius, _delta)
	{
		// can attempt to use array of functors to contain the various shape functions used
		double max_rad_temp = std::max(_radius, d_max);

		//free_mem();

		lx = _lx;
		ly = _ly;
		lz = _lz;

		x_max = _x_inlet + max_rad_temp;
		x_min = _x_inlet - max_rad_temp;
		x_size = std::abs(x_max - x_min);
		y_max = _y_max;
		y_min = _y_min;

		//vol = std::abs((x_max - x_min) * (y_max - y_min) * (z_max - z_min));
		//vol_sqrt = pow(vol, 0.5);

		//std::cout << (z_max - z_min) * (y_max - y_min) << ", " << z_max << ", " << z_min << ", "
		//	<< y_max << ", " << y_min << ", " << ly << ", " << lz << std::endl;
		double yd = y_max - y_min;
		size_t N = trunc((z_max - z_min) * (y_max - y_min) / (4 * ly * lz));
		//size_t N = trunc(std::max((z_max - z_min) / (2 * lz), 1.0) * std::max((y_max - y_min) / (2 * ly), 1.0));
		eddies.resize({ N });
		vf_scaling_factor = 1 / sqrt(N);
		vf_scaling_factor = 1 / pow(N, 0.065);
		//vf_scaling_factor = 1 / sqrt(sqrt((std::max((2*lz)/(z_max-z_min), 1.0) * std::max((2*ly)/(y_max-y_min), 1.0))));
		//vf_scaling_factor = 0.5 * yd < ly ? 2 * sqrt(yd / ly) : 1;
		
		vf_scaling_factor = 1;
		//vf_scaling_factor = 1; // vf scaling factor deals with probability in y, z direction?

		t = 0;

		lt = lx / (2 * _u0); // can use representative if necessary in random distribution?
		t_rand = std::uniform_real_distribution<double>(0, lt);

		instantiate_eddies();
	}

	void increment_eddies(Array<double>& uprime, Array<double>& vprime, Array<double>& wprime,
		const Array<double>& _y_inlet, const Array<double>& _z_inlet, const Array<double>& _a11, const Array<double>& _a21,
		const Array<double>& _a22, const Array<double>& _a31, const Array<double>& _a32, const Array<double>& _a33) {
		t += dt;

		//u_prime.set(0);
		//v_prime.set(0);
		//w_prime.set(0);

		for (size_t idx{ 0 }; idx < eddies.size; idx++) {
			increment_eddy(eddies(idx), uprime, vprime, wprime, _y_inlet, _z_inlet, _a11, _a21, _a22, _a31, _a32, _a33);
			//std::cout << eddies(idx).num_nodes << std::endl;
			//std::cout << "in subregion " << idx << std::endl;
		}
		/*
		for (size_t idx{ 0 }; idx < uprime.size; idx++) {
			uprime(idx) *= vf_scaling_factor;
			vprime(idx) *= vf_scaling_factor;
			wprime(idx) *= vf_scaling_factor;
		}*/
	}

	void increment_eddy(MRSEM_eddy& _eddy, Array<double>& uprime, Array<double>& vprime, Array<double>& wprime,
		const Array<double>& _y_inlet, const Array<double>& _z_inlet, const Array<double>& _a11, const Array<double>& _a21,
		const Array<double>& _a22, const Array<double>& _a31, const Array<double>& _a32, const Array<double>& _a33) {
		bool in_box = _eddy.convect(t);

		if (not in_box) {
			//reset(_eddy);
			reset(_eddy, _y_inlet, _z_inlet);
		}

		for (size_t idx{ 0 }; idx < _eddy.num_nodes; idx++) {
			t_abs = (t - (_eddy.position[0] + _eddy.radius[0])) / _eddy.radius[0];
			y_abs = (_eddy.nodes_pos(idx, 0) - _eddy.position[1]) / _eddy.radius[1];
			z_abs = (_eddy.nodes_pos(idx, 1) - _eddy.position[2]) / _eddy.radius[2];
			r_abs = pow(t_abs * t_abs + y_abs * y_abs + z_abs * z_abs, 0.5);

			if (r_abs < 1) {
				_eddy.shape[0] = shape_t[0](t, _eddy.position[0] + _eddy.radius[0], _eddy.radius[0]);
				_eddy.shape[0] *= shape_t[1](_eddy.nodes_pos(idx, 0), _eddy.position[1], _eddy.radius[1]);
				_eddy.shape[0] *= shape_t[2](_eddy.nodes_pos(idx, 1), _eddy.position[2], _eddy.radius[2]);

				_eddy.shape[1] = shape_y[0](t, _eddy.position[0] + _eddy.radius[0], _eddy.radius[0]);
				_eddy.shape[1] *= shape_y[1](_eddy.nodes_pos(idx, 0), _eddy.position[1], _eddy.radius[1]);
				_eddy.shape[1] *= shape_y[2](_eddy.nodes_pos(idx, 1), _eddy.position[2], _eddy.radius[2]);

				_eddy.shape[2] = shape_z[0](t, _eddy.position[0] + _eddy.radius[0], _eddy.radius[0]);
				_eddy.shape[2] *= shape_z[1](_eddy.nodes_pos(idx, 0), _eddy.position[1], _eddy.radius[1]);
				_eddy.shape[2] *= shape_z[2](_eddy.nodes_pos(idx, 1), _eddy.position[2], _eddy.radius[2]);
			}
			else {
				_eddy.shape[0] = 0;
				_eddy.shape[1] = 0;
				_eddy.shape[2] = 0;
			}


			//uprime(_eddy.nodes(idx, 0), _eddy.nodes(idx, 1)) += vf_scaling_factor * (double)(_eddy.epsilon[0]) * _eddy.shape[0];
			//vprime(_eddy.nodes(idx, 0), _eddy.nodes(idx, 1)) += vf_scaling_factor * (double)(_eddy.epsilon[1]) * _eddy.shape[1];
			//wprime(_eddy.nodes(idx, 0), _eddy.nodes(idx, 1)) += vf_scaling_factor * (double)(_eddy.epsilon[2]) * _eddy.shape[2];

			// edit to include Cholesky decomp terms
			uprime(_eddy.nodes(idx, 0), _eddy.nodes(idx, 1)) += vf_scaling_factor * _a11(_eddy.nodes(idx, 0), _eddy.nodes(idx, 1)) * (double)(_eddy.epsilon[0]) * _eddy.shape[0];
			vprime(_eddy.nodes(idx, 0), _eddy.nodes(idx, 1)) += vf_scaling_factor * (_a21(_eddy.nodes(idx, 0), _eddy.nodes(idx, 1)) * _eddy.epsilon[0] 
				+ _a22(_eddy.nodes(idx, 0), _eddy.nodes(idx, 1)) * (_eddy.epsilon[1])) * _eddy.shape[1];
			wprime(_eddy.nodes(idx, 0), _eddy.nodes(idx, 1)) += vf_scaling_factor * (_a31(_eddy.nodes(idx, 0), _eddy.nodes(idx, 1)) * _eddy.epsilon[0] + 
				_a32(_eddy.nodes(idx, 0), _eddy.nodes(idx, 1)) * _eddy.epsilon[1] + _a33(_eddy.nodes(idx, 0), _eddy.nodes(idx, 1)) * (_eddy.epsilon[2])) * _eddy.shape[2];


			//for (size_t idx2{ 0 }; idx2 < 3; idx2++) {
			//	if (std::abs(_eddy.shape[idx2]) > 1) {
			//		std::cout << t << ", " << idx2 << ", " << _eddy.shape[idx2] << std::endl;
			//	}
			//}
			//std::cout << _eddy.shape[0] << ", " << _eddy.shape[1] << ", " << _eddy.shape[2] << std::endl;
		}
		/*
		for (size_t idx{ 0 }; idx < _y_inlet.size; idx++) {
			// eval normalised variables  (x, y, z) or (t, y, z)
			// eval shape function

			_eddy.shape[0] = shape_t[0](t, _eddy.position[0] + _eddy.radius[0], _eddy.radius[0]);
			_eddy.shape[0] *= shape_t[1](_y_inlet(idx), _eddy.position[1], _eddy.radius[1]);
			_eddy.shape[0] *= shape_t[2](_z_inlet(idx), _eddy.position[2], _eddy.radius[2]);

			_eddy.shape[1] = shape_y[0](t, _eddy.position[0] + _eddy.radius[0], _eddy.radius[0]);
			_eddy.shape[1] *= shape_y[1](_y_inlet(idx), _eddy.position[1], _eddy.radius[1]);
			_eddy.shape[1] *= shape_y[2](_z_inlet(idx), _eddy.position[2], _eddy.radius[2]);

			_eddy.shape[2] = shape_z[0](t, _eddy.position[0] + _eddy.radius[0], _eddy.radius[0]);
			_eddy.shape[2] *= shape_z[1](_y_inlet(idx), _eddy.position[1], _eddy.radius[1]);
			_eddy.shape[2] *= shape_z[2](_z_inlet(idx), _eddy.position[2], _eddy.radius[2]);

			//if (idx == 50) {
				//double tnorm = (t - (_eddy.position[0] + _eddy.radius[0])) / _eddy.radius[0];
				//std::cout << _eddy.shape[0] << ", " << _eddy.shape[1] << ", " << _eddy.shape[2] << ", " << t << ", " << _eddy.position[0] << ", " << _eddy.radius[0] << ", " << tnorm << std::endl;
			//}

			// influence inlet plane, must eval all nodes due to shape function?
			//u_prime(idx) += a11(idx) * _eddy.epsilon(0) * _eddy.shape[0];
			//v_prime(idx) += (a21(idx) * _eddy.epsilon(0)+ a22(idx) * _eddy.epsilon(1)) * _eddy.shape[1];
			//w_prime(idx) += (a31(idx) * _eddy.epsilon(0) + a32(idx) * _eddy.epsilon(1) + a33(idx) * _eddy.epsilon(2))* _eddy.shape[2];

			uprime(idx) += (double)(_eddy.epsilon[0]) * _eddy.shape[0];
			vprime(idx) += (double)(_eddy.epsilon[1]) * _eddy.shape[1];
			wprime(idx) += (double)(_eddy.epsilon[2]) * _eddy.shape[2];
		}*/
	}

private:
	double t_abs;
	double y_abs;
	double z_abs;
	double r_abs;
	double t_temp;
	double y_temp;
	double z_temp;
	double r_temp[3];
	double eps_temp[3];
	double u_temp;

	void instantiate_eddies() {
		for (size_t idx{ 0 }; idx < eddies.size; idx++) {
			reset(eddies(idx));
		}
	}

	const double velocity_fn(const double& y, const double& z) {
		return u0;
	}

	void compute_radius(const double& y, const double& z, const double& u) {
		r_temp[0] = lx / (2 * u0);
		r_temp[1] = ly;
		r_temp[2] = lz;
	}

	void reset(MRSEM_eddy& _eddy) {
		t_temp = t + t_rand(mt);
		//std::cout << "t_temp: " << t << ", " << t_temp << std::endl;
		y_temp = y_rand(mt);
		z_temp = z_rand(mt);
		u_temp = velocity_fn(y_temp, z_temp);
		compute_radius(y_temp, z_temp, u_temp);
		for (size_t e_idx{ 0 }; e_idx < 3; e_idx++) {
			eps_temp[e_idx] = eps_map[eps_rand(mt)];
		}

		_eddy.reset(t_temp, y_temp, z_temp, r_temp, eps_temp, u_temp, dt);//, y_inlet, z_inlet);
	}

	void reset(MRSEM_eddy& _eddy, const Array<double>& _y_inlet, const Array<double>& _z_inlet) {
		t_temp = t + t_rand(mt);
		//std::cout << "t_temp: " << t << ", " << t_temp << std::endl;
		y_temp = y_rand(mt);
		z_temp = z_rand(mt);
		u_temp = velocity_fn(y_temp, z_temp);
		compute_radius(y_temp, z_temp, u_temp);
		for (size_t e_idx{ 0 }; e_idx < 3; e_idx++) {
			eps_temp[e_idx] = eps_map[eps_rand(mt)];
		}

		_eddy.reset(t_temp, y_temp, z_temp, r_temp, eps_temp, u_temp, dt, _y_inlet, _z_inlet);
	}

	void free_mem() {
		x_inlet.resize({ 1 });
		y_inlet.resize({ 1 });
		z_inlet.resize({ 1 });
		a11.resize({ 1 });
		a21.resize({ 1 });
		a22.resize({ 1 });
		a31.resize({ 1 });
		a32.resize({ 1 });
		a33.resize({ 1 });
		u_prime.resize({ 1 });
		v_prime.resize({ 1 });
		w_prime.resize({ 1 });
	}

public:

};
