#pragma once

#include "MRSEM_eddy.h"
#include "region.h"
#include "Array.h"
#include "interpolate.h"
#include "MRSEM_subregion.h"

#include <filesystem>
#include <cstdlib>

// if hairpin vorticies are actual hairpin geometry -> need to create different class for this region specifically?



//-------------------------------------------------------------------

class MRSEM_region {
public:
	Array<MRSEM_subregion> regions;
	Array<double> y_inlet;
	Array<double> z_inlet;
	Array<double> a11;
	Array<double> a21;
	Array<double> a22;
	Array<double> a31;
	Array<double> a32;
	Array<double> a33;

	Array<double> uprime;
	Array<double> vprime;
	Array<double> wprime;

	double u0;

	std::string current_dir = "C:/Users/hayde/Desktop/University_of_Southampton_files/Year_3_(Course_related)/FEEG3003_Individual_Project/Synthetic_eddy_region_code/Synthetic_eddy_cpp/Synthetic_eddy_cpp";

	std::ofstream file;


	MRSEM_region(const Array<double>& _y, const Array<double>& _z, const double& _u0, const size_t& nregions) {
		resize_arrays(_y.size, _z.size);

		Array<double> y_temp;
		Array<double> z_temp;

		y_temp.resize({ _y.size, _z.size });
		z_temp.resize({ _y.size, _z.size });

		y_temp.set(0);
		z_temp.set(0);

		regions.resize({ nregions });
		
		// temporary values
		/*
		for (size_t idx{ 0 }; idx < nregions; idx++) {
			regions(idx) = MRSEM_subregion(1, 1, 1, y_temp, z_temp, 1, 1, 1, 1, 1);
		}*/

		u0 = _u0;

		for (size_t j{ 0 }; j < _y.size; j++) {
			for (size_t k{ 0 }; k < _z.size; k++) {
				y_inlet(j, k) = _y(j);
				z_inlet(j, k) = _z(k);
			}
		}
	}

	void increment_eddies() {
		uprime.set(0);
		vprime.set(0);
		wprime.set(0);

		for (size_t i{ 0 }; i < regions.size; i++) {
			regions(i).increment_eddies(uprime, vprime, wprime, y_inlet, z_inlet);
			//std::cout << uprime(50, 50) << std::endl;
		}
	}

	void set_HIT_RST(const double& TI) {
		for (size_t idx{ 0 }; idx < a11.size; idx++) {
			a11(idx) = u0 * TI;
			a21(idx) = 0;
			a22(idx) = u0 * TI;
			a31(idx) = 0;
			a32(idx) = 0;
			a33(idx) = u0 * TI;
		}
	}
	
	void print_flucts(const size_t& step, const std::string& filename) {
		//std::filesystem::current_path(std::filesystem::temp_directory_path());

		//mkdir

		file.open(current_dir + "/uprime/" + filename + "_uprime_" + std::to_string(step) + ".txt");
		for (size_t j{ 0 }; j < y_inlet.shape[0]; j++) {
			for (size_t k{ 0 }; k < y_inlet.shape[1] - 1; k++) {
				file << uprime(j, k) << ",";
			}
			file << uprime(j, y_inlet.shape[1] - 1) << std::endl;
		}
		file.close();

		file.open(current_dir + "/vprime/" + filename + "_vprime_" + std::to_string(step) + ".txt");
		for (size_t j{ 0 }; j < y_inlet.shape[0]; j++) {
			for (size_t k{ 0 }; k < y_inlet.shape[1] - 1; k++) {
				file << vprime(j, k) << ",";
			}
			file << vprime(j, y_inlet.shape[1] - 1) << std::endl;
		}
		file.close();

		file.open(current_dir + "/wprime/" + filename + "_wprime_" + std::to_string(step) + ".txt");
		for (size_t j{ 0 }; j < y_inlet.shape[0]; j++) {
			for (size_t k{ 0 }; k < y_inlet.shape[1] - 1; k++) {
				file << wprime(j, k) << ",";
			}
			file << wprime(j, y_inlet.shape[1] - 1) << std::endl;
		}
		file.close();
	}
	

private:
	void resize_arrays(const size_t& ysize, const size_t& zsize) {
		y_inlet.resize({ ysize, zsize });
		z_inlet.resize({ ysize, zsize });

		a11.resize({ ysize, zsize });
		a21.resize({ ysize, zsize });
		a22.resize({ ysize, zsize });
		a31.resize({ ysize, zsize });
		a32.resize({ ysize, zsize });
		a33.resize({ ysize, zsize });

		uprime.resize({ ysize, zsize });
		vprime.resize({ ysize, zsize });
		wprime.resize({ ysize, zsize });
	}
};

