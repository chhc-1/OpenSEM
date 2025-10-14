#pragma once

#include "Array.h"
#include "Eddy.h"
#include "region.h"
#include "DFSEM_eddy.h"
#include <cstdlib>
#include <random>


class DFSEM_region : public region {
public:
	Array<DFSEM_eddy> eddies;
	double C1;
	double C2;
	double avg_radius;

	// eigenvectors and values for each grid point
	Array<double> eigval1;
	Array<double> eigval2;
	Array<double> eigval3;

	Array<double> eigvect11;
	Array<double> eigvect12;
	Array<double> eigvect13;
	Array<double> eigvect21;
	Array<double> eigvect22;
	Array<double> eigvect23;
	Array<double> eigvect31;
	Array<double> eigvect32;
	Array<double> eigvect33;

	/*
	Array<double> inv_eigvect11;
	Array<double> inv_eigvect12;
	Array<double> inv_eigvect13;
	Array<double> inv_eigvect21;
	Array<double> inv_eigvect22;
	Array<double> inv_eigvect23;
	Array<double> inv_eigvect31;
	Array<double> inv_eigvect32;
	Array<double> inv_eigvect33;
	*/
	Array<double> u_prime_local;
	Array<double> v_prime_local;
	Array<double> w_prime_local;

	// constructor here
	
	DFSEM_region(const double& _u0, const double& _dt, const double& _x_inlet, const Array<double>& _y_inlet,
		const Array<double>& _z_inlet, const double& _min_radius, const double& _rep_radius, const double& _max_radius, const double& _delta)
		: region(_u0, _dt, _x_inlet, _y_inlet, _z_inlet, _rep_radius, _delta)
	{
		x_max = _x_inlet + _max_radius;
		x_min = _x_inlet - _max_radius;
		x_size = x_max - x_min;

		vol = std::abs((x_max - x_min) * (y_max - y_min) * (z_max - z_min));
		vol_sqrt = pow(vol, 0.5);

		size_t N = trunc(vol / pow(_rep_radius, 3));
		eddies.resize({ N });

		vf_scaling_factor = 1 / pow(N, 0.5);
		avg_radius = _rep_radius;

		RST_idx[0] = 0;
		RST_idx[1] = 2;
		RST_idx[2] = 5;

		calc_C1(_min_radius);

		calc_C2();

		resize_eigs();

		instantiate_eigs();

		instantiate_eddies();
	}


	void increment_eddies() {
		u_prime_local.set(0);
		v_prime_local.set(0);
		w_prime_local.set(0);

		for (size_t idx{ 0 }; idx < eddies.size; idx++) {
			increment_eddy(eddies(idx));
		}

		for (size_t i{ 0 }; i < u_prime.size; i++) {
			u_prime_local(i) *= vf_scaling_factor;
			v_prime_local(i) *= vf_scaling_factor;
			w_prime_local(i) *= vf_scaling_factor;

			//std::cout << u_prime_local(i) << ", ";
		}
		
		loc_to_gbl();
	}

	void increment_eddy(DFSEM_eddy& _eddy) {
		in_box = _eddy.convect(x_max);
		if (not in_box) {
			y_temp = y_rand(mt);
			z_temp = z_rand(mt);
			radius_fn(y_temp, z_temp);
			RST_calc(y_temp, z_temp);
			RST_calc_eigval(RST_temp);
			alpha_fn(radius_temp);

			// can change x_min to reduce eddy x position instead
			// also consider removing y_inlet and z_inlet - do not store the nodes to iterate through in each eddy
			// can also add function for eddy velocity
			_eddy.reset(x_min, y_temp, z_temp, radius_temp, vol, u0, dt, y_inlet, z_inlet, alpha_temp);
		}
		double temp_x = x_inlet(0, 0);
		for (size_t i{ 0 }; i < _eddy.num_nodes; i++) {
			_eddy.calc_r(temp_x, _eddy.nodes_pos(i, 0), _eddy.nodes_pos(i, 1));
			u_prime_local(_eddy.nodes(i, 0), _eddy.nodes(i, 1)) += _eddy.radius[0] * (1 - _eddy.r_magn2) * (_eddy.rk[1] * _eddy.alpha[2] - _eddy.rk[2] * _eddy.alpha[1]);
			v_prime_local(_eddy.nodes(i, 0), _eddy.nodes(i, 1)) += _eddy.radius[1] * (1 - _eddy.r_magn2) * (_eddy.rk[2] * _eddy.alpha[0] - _eddy.rk[0] * _eddy.alpha[2]);
			w_prime_local(_eddy.nodes(i, 0), _eddy.nodes(i, 1)) += _eddy.radius[2] * (1 - _eddy.r_magn2) * (_eddy.rk[0] * _eddy.alpha[1] - _eddy.rk[1] * _eddy.alpha[0]);
		}

		/*
		for (size_t i{ 0 }; i < u_prime.size; i++) {
			_eddy.calc_r(x_inlet(i), y_inlet(i), z_inlet(i));

			//std::cout << _eddy.radius[0] << ", " << _eddy.r_magn2 << std::endl;
			//std::cout << _eddy.alpha[0] << std::endl;

			u_prime_local(i) += _eddy.radius[0] * (1 - _eddy.r_magn2) * (_eddy.rk[1] * _eddy.alpha[2] - _eddy.rk[2] * _eddy.alpha[1]);
			v_prime_local(i) += _eddy.radius[1] * (1 - _eddy.r_magn2) * (_eddy.rk[2] * _eddy.alpha[0] - _eddy.rk[0] * _eddy.alpha[2]);
			w_prime_local(i) += _eddy.radius[2] * (1 - _eddy.r_magn2) * (_eddy.rk[0] * _eddy.alpha[1] - _eddy.rk[1] * _eddy.alpha[0]);
		}
		*/
	}

	// transfers local velocity fluctuation to global velocity fluctuation
	void loc_to_gbl() {
		for (size_t j{ 0 }; j < u_prime.shape[0]; j++) {
			for (size_t k{ 0 }; k < u_prime.shape[1]; k++) {
				u_prime(j, k) = C1 * (eigvect11(j, k) * u_prime_local(j, k) + eigvect12(j, k) * v_prime_local(j, k) + eigvect13(j, k) * w_prime_local(j, k));
				v_prime(j, k) = C1 * (eigvect21(j, k) * u_prime_local(j, k) + eigvect22(j, k) * v_prime_local(j, k) + eigvect23(j, k) * w_prime_local(j, k));
				w_prime(j, k) = C1 * (eigvect31(j, k) * u_prime_local(j, k) + eigvect32(j, k) * v_prime_local(j, k) + eigvect33(j, k) * w_prime_local(j, k));
			}
		}
	}

private:
	double RST_temp[6]; // returns uu, uv, vv, uw, vw, ww for arbitrary position
	size_t RST_idx[3];
	double RST_eigval_temp[3];
	double e_matrix[9];
	double RHS_vect[3];
	double eigvect_temp[9]; // returns eigenvectors for given 3x3 matrix (RST), 0-2 is 1st eigvect, 3-5 is 2nd eigvect, 6-8 is 3rd eigvect
	double radius_temp[3];
	double alpha_var[3];
	double alpha_temp[3];
	double sum_temp;
	double a, b, c, d;
	double A, Btilde, C, C3, D;
	double real_root_check;
	double c_val1, c_val2;
	double theta;
	//double alpha_max[3];

	void radius_fn(const double& y, const double& z) {
		// eval radius_temp
		radius_temp[0] = avg_radius;
		radius_temp[1] = avg_radius;
		radius_temp[2] = avg_radius;
	}

	void alpha_fn(const double _radius[3]) {
		sum_temp = 0;
		for (size_t c{ 0 }; c < 3; c++) {
			sum_temp += RST_temp[RST_idx[c]] / (_radius[c] * _radius[c]);
		}

		// need to replace 1e-7 with a different truncation value
		for (size_t c{ 0 }; c < 3; c++) {
			alpha_var[c] = std::max(1e-7, (sum_temp - 2 * RST_temp[RST_idx[c]] / (_radius[c] * _radius[c])) / (2 * C2));
			// line below may not be correct...
			alpha_temp[c] = alpha_var[c] * double(eps_map[eps_rand(mt)]);
			
		}
	}

	double* RST_calc(const double& y, const double& z) {
		// RST_temp[0] = 
		// RST_temp[1] =
		
		// temporary HIT RST for 1% turbulence intensity
		RST_temp[0] = 0.01 * u0 * 0.01 * u0;
		RST_temp[1] = 0;
		RST_temp[2] = 0.01 * u0 * 0.01 * u0;
		RST_temp[3] = 0;
		RST_temp[4] = 0;
		RST_temp[5] = 0.01 * u0 * 0.01 * u0;
		return &RST_temp[0];
	}

	void RST_calc_eigval(const double RST[6]) {
		// cubic eval for eigenvalues for Reynold's Stress Tensor
		// may need to give warning for imaginary eigenvalues for RST?
		a = 1;
		b = -(RST[0] + RST[2] + RST[5]);
		c = RST[0] * RST[2] + RST[0] * RST[5] + RST[2] * RST[5] + RST[1] * RST[1] + RST[3] * RST[3] + RST[4] * RST[4];
		d = -(RST[0] * RST[2] * RST[5] - RST[0] * RST[4] * RST[4] - RST[1] * RST[1] * RST[5]
			+ RST[1] * RST[3] * RST[4] + RST[4] * RST[1] + RST[3] - RST[2] * RST[3] * RST[3]);

		//c_sqrt_val1 = (2 * b * b * b - 9 * a * b * c + 27 * a * a * d);
		real_root_check = -27 * a * a * d * d + 18 * a * b * c * d - 4 * a * c * c * c - 4 * b * b * b * d + b * b * c * c;
		if (real_root_check < 0) {
			std::cout << "WARNING: c_sqrt_val2 less than 0. Cannot obtain real eigvals. proceed?" << std::endl;
			std::cin.get();
		}
		A = 2 * b * b * b - 9 * a * b * c + 27 * a * a * d;
		Btilde = pow(4 * pow(b * b - 3 * a * c, 3) - pow(2 * b * b * b - 9 * a * b * c + 27 * a * a * d, 2), 0.5);
		D = -b / (3*a);
		C = pow(A * A + Btilde * Btilde, 0.5);
		C3 = cbrt(C/2);
		theta = atan2(Btilde, A);
		RST_eigval_temp[0] = D - 2 / (3 * a) * C3 * cos(theta / 3);
		RST_eigval_temp[1] = D + 2 / (3 * a) * C3 * cos(theta / 3 + M_PI / 3);
		RST_eigval_temp[2] = D + 2 / (3 * a) * C3 * cos(theta / 3 - M_PI / 3);
		/*
		for (size_t o{ 0 }; o < 3; o++) {
			std::cout << RST_eigval_temp[o] << ", ";
		}*/
	}

	void RST_calc_eigvect(const double RST[6], const double eigval[3]) {
		/*
		for (size_t i{ 0 }; i < 3; i++) {
			e_matrix[0] = RST[0] - eigval[i];
			e_matrix[1] = RST[1];
			e_matrix[2] = RST[3];
			e_matrix[3] = RST[1];
			e_matrix[4] = RST[2] - eigval[i];
			e_matrix[5] = RST[4];
			e_matrix[6] = RST[3];
			e_matrix[7] = RST[4];
			e_matrix[8] = RST[5] - eigval[i];

			RHS_vect[0] = 0;

			// for HIT, RST is diagonal - e_matrix is all 0...
			for (size_t l{ 1 }; l < 3; l++) {
				for (size_t m{ l }; m < 3; m++) {
					for (size_t n{ m - 1 }; n < 3; n++) {
						e_matrix[3 * m + n] -= e_matrix[3 * m + m - 1] / e_matrix[3 * (l - 1) + l - 1] * e_matrix[3 * (l - 1) + l - 1 + n];
					}
					RHS_vect[m] -= RHS_vect[l - 1] * e_matrix[3 * m + m - 1] / e_matrix[3 * (l - 1) + l - 1];
				}
			}

			for (size_t m{ 0 }; m < 3; m++) {
				for (size_t n{ 0 }; n < 3; n++) {
					std::cout << e_matrix[3 * m + n] << ", ";
				}
				std::cout << std::endl;
			}

			for (size_t m{ 2 }; m > 0; m--) {
				eigvect_temp[3 * i + m] = RHS_vect[m];
				for (size_t n{ 2-m }; n < 2; n--) {
					eigvect_temp[3 * i + m] -= e_matrix[3 * m + n] * eigvect_temp[3 * i + m + n];
				}
				eigvect_temp[3 * i + m] /= e_matrix[3 * m + m];
			}
		}

		for (size_t m{ 0 }; m < 3; m++) {
			for (size_t n{ 0 }; n < 3; n++) {
				std::cout << eigvect_temp[3 * m + n] << ", ";
			}
			std::cout << std::endl;
		}
		*/
		for (size_t i{ 0 }; i < 3; i++) {
			for (size_t m{ 0 }; m < 3; m++) {
				if (i == m) {
					eigvect_temp[3 * i + m] = 1;
				}
				else {
					eigvect_temp[3 * i + m] = 0;
				}
			}
		}

	}

	void instantiate_eddies() {
		std::uniform_real_distribution<double> x_rand = std::uniform_real_distribution<double>(x_min, x_max);
			
		for (size_t i{ 0 }; i < eddies.size; i++) {
			eddies(i) = DFSEM_eddy(100); // 0
			x_temp = x_rand(mt);
			y_temp = y_rand(mt);
			z_temp = z_rand(mt);

			radius_fn(y_temp, z_temp);
			RST_calc(y_temp, z_temp);
			RST_calc_eigval(RST_temp);
			alpha_fn(radius_temp);

			eddies(i).reset(x_temp, y_temp, z_temp, radius_temp, vol, u0, dt, y_inlet, z_inlet, alpha_temp);
		}
	}
	void resize_eigs() {
		eigval1.resize({ y_inlet.shape[0], y_inlet.shape[1] });
		eigval2.resize({ y_inlet.shape[0], y_inlet.shape[1] });
		eigval3.resize({ y_inlet.shape[0], y_inlet.shape[1] });

		eigvect11.resize({ y_inlet.shape[0], y_inlet.shape[1] });
		eigvect12.resize({ y_inlet.shape[0], y_inlet.shape[1] });
		eigvect13.resize({ y_inlet.shape[0], y_inlet.shape[1] });
		eigvect21.resize({ y_inlet.shape[0], y_inlet.shape[1] });
		eigvect22.resize({ y_inlet.shape[0], y_inlet.shape[1] });
		eigvect23.resize({ y_inlet.shape[0], y_inlet.shape[1] });
		eigvect31.resize({ y_inlet.shape[0], y_inlet.shape[1] });
		eigvect32.resize({ y_inlet.shape[0], y_inlet.shape[1] });
		eigvect33.resize({ y_inlet.shape[0], y_inlet.shape[1] });

		u_prime_local.resize({ y_inlet.shape[0], y_inlet.shape[1] });
		v_prime_local.resize({ y_inlet.shape[0], y_inlet.shape[1] });
		w_prime_local.resize({ y_inlet.shape[0], y_inlet.shape[1] });
	}

	void instantiate_eigs() {
		for (size_t j{ 0 }; j < y_inlet.shape[0]; j++) {
			for (size_t k{ 0 }; k < y_inlet.shape[1]; k++) {
				RST_calc(y_inlet(j, k), z_inlet(j, k));
				RST_calc_eigval(RST_temp);

				eigval1(j, k) = RST_eigval_temp[0];
				eigval2(j, k) = RST_eigval_temp[1];
				eigval3(j, k) = RST_eigval_temp[2];

				RST_calc_eigvect(RST_temp, RST_eigval_temp);

				eigvect11(j, k) = eigvect_temp[0];
				eigvect12(j, k) = eigvect_temp[1];
				eigvect13(j, k) = eigvect_temp[2];
				eigvect21(j, k) = eigvect_temp[3];
				eigvect22(j, k) = eigvect_temp[4];
				eigvect23(j, k) = eigvect_temp[5];
				eigvect31(j, k) = eigvect_temp[6];
				eigvect32(j, k) = eigvect_temp[7];
				eigvect33(j, k) = eigvect_temp[8];
			}
		}
	}

	void calc_C1(const double& min_rad) {
		C1 = sqrt(10 * vol) * (avg_radius) / (sqrt(eddies.size) * (pow(avg_radius, 3))) * avg_radius;
	}

	void calc_C2() {
		// C2 scaling factor is currently INCORRECT - needs fix
		C2 = 2.0 / sqrt(3.6);
	}


public:


};


/*
class DFSEM_region : public region{
public:
	Array<DFSEM_eddy> eddies;
	double min_radius;
	std::normal_distribution<double> radius_dist;
	double C2; // need method to calculate C2
	double C1; // use representative radius for calculating C1

	Array<double> eigval1;
	Array<double> eigval2;
	Array<double> eigval3;

	Array<double> eigvect11;
	Array<double> eigvect12;
	Array<double> eigvect13;
	Array<double> eigvect21;
	Array<double> eigvect22;
	Array<double> eigvect23;
	Array<double> eigvect31;
	Array<double> eigvect32;
	Array<double> eigvect33;

	// can use random distribution for norm of 0 - 1 numbers, then multiply by scaling to achieve new values for radius/radii


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

			_eddy.reset(_eddy.position(0) - x_size, y_temp, z_temp, rad_vect_temp.array1, vol, u0, dt, y_inlet, z_inlet, eps_temp);
		}

		double temp_x = x_inlet(0, 0); // assume inlet is planar constant x


		for (size_t i{ 0 }; i < _eddy.num_nodes; i++) {
			_eddy.calc_r(temp_x, _eddy.nodes_pos(i, 0), _eddy.nodes_pos(i, 1));
			//_eddy.shape = _eddy.shape_scaling_factor;
			//_eddy.shape *= _eddy.shape_fn(temp_x, _eddy.nodes_pos(i, 0), _eddy.nodes_pos(i, 1));

			//u_prime(_eddy.nodes(i, 0), _eddy.nodes(i, 1)) += (_eddy.rk(1) * _eddy.epsilon(2) - _eddy.rk(2) * _eddy.epsilon(1)) * _eddy.shape / _eddy.r_magn3;
			//v_prime(_eddy.nodes(i, 0), _eddy.nodes(i, 1)) += (_eddy.rk(2) * _eddy.epsilon(0) - _eddy.rk(0) * _eddy.epsilon(2)) * _eddy.shape / _eddy.r_magn3;
			//w_prime(_eddy.nodes(i, 0), _eddy.nodes(i, 1)) += (_eddy.rk(0) * _eddy.epsilon(1) - _eddy.rk(1) * _eddy.epsilon(0)) * _eddy.shape / _eddy.r_magn3;

			u_prime(_eddy.nodes(i, 0), _eddy.nodes(i, 1)) += (_eddy.rk[1] * _eddy.epsilon(2) - _eddy.rk[2] * _eddy.epsilon(1)) * _eddy.shape_fn(_eddy.radius[0]);
			v_prime(_eddy.nodes(i, 0), _eddy.nodes(i, 1)) += (_eddy.rk[2] * _eddy.epsilon(0) - _eddy.rk[0] * _eddy.epsilon(2)) * _eddy.shape_fn(_eddy.radius[1]);
			w_prime(_eddy.nodes(i, 0), _eddy.nodes(i, 1)) += (_eddy.rk[0] * _eddy.epsilon(1) - _eddy.rk[1] * _eddy.epsilon(0)) * _eddy.shape_fn(_eddy.radius[2]);
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

			eddies(i).reset(x_temp, y_temp, z_temp, rad_vect_temp.array1, vol, u0, dt, y_inlet, z_inlet, eps_temp);

			eddies(i).K = 0.01;//1e-5; /// temporary - needs normalisation
		}
	}

	double mean_temp_radius;
	double radius_std;
	//double radius_temp;
	Array<double> rad_vect_temp;

public:



};
*/


