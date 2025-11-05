#pragma once

#include "eddy.h"
#include "Array.h"
#include <cstdlib>

#define M_PI 3.14159265358979

struct MRSEM_eddy : public eddy {

	double shape[3];
	double radius[3]; // first axis is time, not streamwise position
	double epsilon[3];
	double position[3];
	// position(0) is now time position when eddy is created

	MRSEM_eddy()
		:eddy(300)
	{

	}

	const bool convect(const double& t) {
		if ((t - position[0] - radius[0]) / radius[0] > 1) {
			return false;
		}
		else {
			return true;
		}
	}

	void reset(const double& x_new, const double& y_new, const double& z_new, double* r_new, double* eps_new,
		const double& u, const double& dt){//, const Array<double>& y_plane, const Array<double>& z_plane) {
		position[0] = x_new;
		position[1] = y_new;
		position[2] = z_new;

		increment = u * dt;
		num_nodes = 0;

		for (size_t e_idx{ 0 }; e_idx < 3; e_idx++) {
			radius[e_idx] = r_new[e_idx];
			epsilon[e_idx] = eps_new[e_idx];
		}
		
		double dy;
		double dz;

		/*
		for (size_t j{ 0 }; j < y_plane.shape[0]; j++) {
			for (size_t k{ 0 }; k < y_plane.shape[1]; k++) {

				// eval delta y, delta z instead of evaluating radius?
				dy = y_plane(j, k) - radius[1];
				dz = z_plane(j, k) - radius[2];

				if (std::abs(dy) < radius[1] && std::abs(dz) < radius[2]) {
					nodes(num_nodes, 0) = j;
					nodes(num_nodes, 1) = k;
					nodes_pos(num_nodes, 0) = y_plane(j, k);
					nodes_pos(num_nodes, 1) = z_plane(j, k);
					num_nodes++;
				}
			}
		}*/
	}

	void reset(const double& x_new, const double& y_new, const double& z_new, double* r_new, double* eps_new,
		const double& u, const double& dt, const Array<double>& y_plane, const Array<double>& z_plane) {
		position[0] = x_new;
		position[1] = y_new;
		position[2] = z_new;

		increment = u * dt;
		num_nodes = 0;

		for (size_t e_idx{ 0 }; e_idx < 3; e_idx++) {
			radius[e_idx] = r_new[e_idx];
			epsilon[e_idx] = eps_new[e_idx];
		}

		double dy;
		double dz;

		
		for (size_t j{ 0 }; j < y_plane.shape[0]; j++) {
			for (size_t k{ 0 }; k < y_plane.shape[1]; k++) {

				// eval delta y, delta z instead of evaluating radius?
				dy = y_plane(j, k) - position[1];
				dz = z_plane(j, k) - position[2];

				if (std::abs(dy) < radius[1] && std::abs(dz) < radius[2]) {
					nodes(num_nodes, 0) = j;
					nodes(num_nodes, 1) = k;
					nodes_pos(num_nodes, 0) = y_plane(j, k);
					nodes_pos(num_nodes, 1) = z_plane(j, k);
					num_nodes++;
				}
			}
		}
	}

};

struct shape_fn {
	double (*f)(const double&) = nullptr; // store function pointer

	shape_fn() {

	}

	shape_fn(double (*_f)(const double&)) {
		f = _f;
	}

	virtual double eval(const double& xtilde) {
		if (std::abs(xtilde) < 1) {
			return f(xtilde);
		}
		else {
			return 0;
		}
	}
	virtual double eval(const double& x_in, const double& x_self, const double& radius) {
		double x_norm = (x_in - x_self) / radius;
		if (std::abs(x_norm) < 1) {
			return f(x_norm);
		}
		else {
			return 0;
		}
	}

	virtual double operator()(const double& xtilde) {
		return eval(xtilde);
	}
	virtual double operator()(const double& x_in, const double& x_self, const double& radius) {
		return eval(x_in, x_self, radius);
	}
};

/*
struct gaussian_fn : public shape_fn {

	// use unit variance, zero mean function
	double eval(const double& xtilde) override {
		if (std::abs(xtilde) < 1) {
			return 1 / sqrt(2 * M_PI) * exp(-0.5 * xtilde);
		}
		else {
			return 0;
		}
	}
	double eval(const double& x_in, const double& x_self, const double& radius) override {
		double x_norm = (x_in - x_self) / radius;
		//std::cout << x_norm << " gaussian" << std::endl;
		if (std::abs(x_norm) < 1) {
			return 1 / sqrt(2 * M_PI) * exp(-0.5 * x_norm);
		}
		else {
			return 0;
		}
	}

	double operator()(const double& xtilde) override {
		return eval(xtilde);
	}
	double operator()(const double& x_in, const double& x_self, const double& radius) override  {
		return eval(x_in, x_self, radius);
	}
};

struct mgaussian_fn : public shape_fn {

	// use unit variance, zero mean function
	double eval(const double& xtilde) {
		if (std::abs(xtilde) < 1) {
			return -1 / sqrt(2 * M_PI) * exp(-0.5 * xtilde);
		}
		else {
			return 0;
		}
	}
	double eval(const double& x_in, const double& x_self, const double& radius) {
		double x_norm = (x_in - x_self) / radius;
		if (std::abs(x_norm) < 1) {
			return -1 / sqrt(2 * M_PI) * exp(-0.5 * (x_in - x_self) / radius);
		}
		else {
			return 0;
		}
	}

	double operator()(const double& xtilde) {
		return eval(xtilde);
	}
	double operator()(const double& x_in, const double& x_self, const double& radius) {
		return eval(x_in, x_self, radius);
	}
};

struct xgaussian_fn : public shape_fn {
	double eval(const double& xtilde) {
		if (std::abs(xtilde) < 1) {
			return xtilde / sqrt(2 * M_PI) * exp(-0.5 * xtilde);
		}
		else {
			return 0;
		}
	}
	double eval(const double& x_in, const double& x_self, const double& radius) {
		double x_norm = (x_in - x_self) / radius;
		if (std::abs(x_norm) < 1) {
			return x_norm / sqrt(2 * M_PI) * exp(-0.5 * x_norm);
		}
		else {
			return 0;
		}
	}

	double operator()(const double& xtilde) {
		return eval(xtilde);
	}
	double operator()(const double& x_in, const double& x_self, const double& radius) {
		return eval(x_in, x_self, radius);
	}
};

struct H_fn : public shape_fn {

	double eval(const double& xtilde) {
		if (std::abs(xtilde) < 1) {
			return 1 - cos(2 * M_PI * xtilde) / (2 * M_PI * xtilde * sqrt(0.214));
		}
		else {
			return 0;
		}
	}
	double eval(const double& x_in, const double& x_self, const double& radius) {
		double x_norm = (x_in - x_self) / radius;
		if (std::abs(x_norm) < 1) {
			return 1 - cos(2 * M_PI * x_norm) / (2 * M_PI * x_norm * sqrt(0.214));
		}
		else {
			return 0;
		}
	}

	double operator()(const double& xtilde) {
		return eval(xtilde);
	}
	double operator()(const double& x_in, const double& x_self, const double& radius) {
		return eval(x_in, x_self, radius);
	}
};*/