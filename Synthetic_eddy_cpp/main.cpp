#include "Array.h"
#include "Eddy.h"
#include "oSEM_eddy.h"
#include "oSEM_region.h"
#include "ISEM1_eddy.h"
#include "ISEM1_region.h"
#include "DFSEM_eddy.h"
#include "DFSEM_region.h"
#include "MRSEM_eddy.h"
#include "MRSEM_subregion.h"
#include "MRSEM_region.h"
#include "TBL_MRSEM.h"
#include "interpolate.h"
#include "target_stats.h"

#include <cstdlib>
#include <random>
#include <map>
#include <string>
#include <iomanip>
#include <chrono>
#include <fstream>

struct AllocationTracker {
public:
	size_t current_allocs{ 0 };
	size_t num_allocs{ 0 };
	size_t deleted_allocs{ 0 };

	void print_allocs() {
		std::cout << "num_allocs: " << num_allocs << " | current_allocs: " << current_allocs << std::endl;
	}

};

static AllocationTracker _tracker;

void* operator new(size_t size) {
	_tracker.num_allocs += 1;
	_tracker.current_allocs += size;
	return malloc(size);
}

/*
void operator delete(void* _ptr) {
	_tracker.deleted_allocs++;
	//_tracker.current_allocs -= ;
}*/
/*
struct timer {
	std::chrono::high_resolution_clock::time_point start_time;
	std::chrono::high_resolution_clock::time_point end_time;
	std::chrono::duration<double> dt;

	timer() {
		start_time = std::chrono::high_resolution_clock::now();
	}
	~timer() {
		end_time = std::chrono::high_resolution_clock::now();
		dt = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
		std::cout << dt.count() << " seconds" << std::endl;
	}
};*/



int main(){

	std::random_device rd;
	std::uniform_real_distribution<double> dist1(1, 100);
	//std::mt19937 mt(rd());
	std::mt19937 mt(1729);

	std::uniform_int_distribution<int> dist2(0, 1);
	std::map<int, int> map1{ {0, -1}, {1, 1}};

	/*
	for (size_t i{ 0 }; i < 100; i++) {
		std::cout << dist1(mt) << std::endl;
	}
	
	std::cout << "-----------------------" << std::endl;

	for (size_t i{ 0 }; i < 100; i++) {
		std::cout << map1[dist2(mt)] << std::endl;
	}
	*/

	double delta = 0.007;

	double x_pos = 0.0;
	double y_min = 0.0;
	double y_max = 1.4 * delta;
	double z_min = 0.0;
	double z_max = 0.05;

	size_t n_y = 100; // 200;
	size_t n_z = 150;

	Array<double> y_pos;
	Array<double> z_pos;

	y_pos.LinRange(y_min, y_max, n_y);
	z_pos.LinRange(z_min, z_max, n_z);

	Array<double> x_plane;
	Array<double> y_plane;
	Array<double> z_plane;

	x_plane.resize({ n_y, n_z });
	y_plane.resize({ n_y, n_z });
	z_plane.resize({ n_y, n_z });

	x_plane.set(x_pos);
	for (size_t j{ 0 }; j < n_y; j++) {
		for (size_t k{ 0 }; k < n_z; k++) {
			y_plane(j, k) = y_pos(j);
			z_plane(j, k) = z_pos(k);
		}
	}

	//double min_rad1 = 0.05 * delta;
	double rep_rad1 = 0.2 * delta; // 0.2 * delta;
	double min_rad1 = 0.1 * delta;
	//double max_rad1 = 0.286 * delta;
	double max_rad1 = 0.41 * delta;

	double u0 = 823.6; // 100;
	double utau = 40.6; // utau = 40.6m/s
	//double u0 = 150;
	double dt = 0.000002;
	double nu = 0.0006129165;


	//std::cout << region.C1 << std::endl; // DFSEM only


	//_tracker.print_allocs();
	//region.increment_eddies();
	//_tracker.print_allocs();
	//region.increment_eddy(region.eddies(0));
	//_tracker.print_allocs();

	//int test1;

	
	//test1 = dist3(urng);
	//_tracker.print_allocs();

	//std::cout << region.a11(0) << std::endl;

	//_tracker.print_allocs();
	//region.increment_eddies();
	//_tracker.print_allocs();
	//Array<double> arr1;
	//arr1.LinRange(0, 2 * M_PI, 1000);
	//Array<double> arr2;
	//arr2.LinRange(0, 4 * M_PI, 1000);

	//interpolator<double> interp1 = interpolator<double>(arr1, arr2);

	//std::cout << interp1(M_PI) << std::endl;

	size_t dpoints{260};
	Array<double> y_interp_arr = Array<double>({ dpoints });
	Array<double> u_interp_arr = Array<double>({ dpoints });
	Array<double> uu_interp_arr = Array<double>({ dpoints });
	Array<double> vv_interp_arr = Array<double>({ dpoints });
	Array<double> ww_interp_arr = Array<double>({ dpoints });
	Array<double> uv_interp_arr = Array<double>({ dpoints });

	std::cout << dpoints << std::endl;

	y_interp_arr.array1 = y_inp;
	u_interp_arr.array1 = uinf_inp;
	uu_interp_arr.array1 = uu_inp;
	vv_interp_arr.array1 = vv_inp;
	ww_interp_arr.array1 = ww_inp;
	uv_interp_arr.array1 = uv_inp;


	/*
	for (size_t i{ 0 }; i < dpoints; i++) {
		y_interp_arr(i) = y_inp[i];
		uu_interp_arr(i) = uu_inp[i];
		vv_interp_arr(i) = vv_inp[i];
		ww_interp_arr(i) = ww_inp[i];
		uv_interp_arr(i) = uv_inp[i];
	}
	*/

	interpolator<double> u_interp(y_interp_arr, u_interp_arr);
	interpolator<double> uu_interp(y_interp_arr, uu_interp_arr);
	interpolator<double> vv_interp(y_interp_arr, vv_interp_arr);
	interpolator<double> ww_interp(y_interp_arr, ww_interp_arr);
	interpolator<double> uv_interp(y_interp_arr, uv_interp_arr);

	Array<interpolator<double>> interps_arr;
	interps_arr.resize({ 4 });
	interps_arr(0) = uu_interp;
	interps_arr(1) = uv_interp;
	interps_arr(2) = vv_interp;
	interps_arr(3) = ww_interp;

	//ISEM1_region region = ISEM1_region(u0, dt, x_pos, y_plane, z_plane, rep_rad1, max_rad1, delta, u_interp);

	//y_min -= max_rad1;
	//oSEM_region region = oSEM_region(u0, dt, x_pos, y_plane, z_plane, max_rad1, delta);

	//DFSEM_region region = DFSEM_region(u0, dt, x_pos, y_plane, z_plane, min_rad1, rep_rad1, max_rad1, delta,
	//	interps_arr);

	dt = 5e-7;
	MRSEM_region region = TBL_MRSEM(n_y, n_z, y_max, z_max, u0, dt, delta, utau, nu);

	//std::cout << region.eddies.size << std::endl;

	double TI1 = 0.01;

	region.set_HIT_RST(TI1);

	//region.set_RST_interps(uu_interp, uv_interp, vv_interp, ww_interp);

	Array<double> r11({ n_y, n_z });
	Array<double> r21({ n_y, n_z });
	Array<double> r22({ n_y, n_z });
	Array<double> r31({ n_y, n_z });
	Array<double> r32({ n_y, n_z });
	Array<double> r33({ n_y, n_z });
	
	for (size_t j{ 0 }; j < n_y; j++) {
		for (size_t k{ 0 }; k < n_z; k++) {
			r11(j, k) = uu_interp(y_pos(j));
			r21(j, k) = uv_interp(y_pos(j));
			r22(j, k) = vv_interp(y_pos(j));
			r31(j, k) = 0;
			r32(j, k) = 0;
			r33(j, k) = ww_interp(y_pos(j));
		}

		//std::cout << y_pos(j) << ", " << uu_interp(y_pos(j)) << std::endl;
	}

	//region.set_RST(r11, r21, r22, r31, r32, r33);
	//region.dt = 2e-8;
	// tke also varies with wall distance -> fixed?
	

	//inp_data.close();

	size_t iters{ 300 };

	std::ofstream file;
	
	file.open("RST.txt");
	for (size_t i{ 0 }; i < n_y; i++) {
		file << r11(i, 0) << ", " << r21(i, 0) << ", " << r22(i, 0) << ", " << r31(i, 0) << ", " << r32(i, 0) << ", " << r33(i, 0) << std::endl;
	}
	file.close();
	/*
	for (size_t idx{ 0 }; idx < n_y; idx++) {
		std::cout << region.a11(idx, 0) << ", " << region.a21(idx, 0) << ", " << region.a22(idx, 0) << ", " << region.a33(idx, 0) << ", " << std::endl;
	}*/

	/*
	for (size_t j{ 0 }; j < region.a11.shape[0]; j++) {
		for (size_t k{ 0 }; k < region.a11.shape[1]; k++) {
			std::cout << region.a11(j, k) << ", ";
			std::cout << region.a21(j, k) << ", ";
			std::cout << region.a22(j, k) << ", ";
			std::cout << region.a31(j, k) << ", ";
			std::cout << region.a32(j, k) << ", ";
			std::cout << region.a33(j, k) << std::endl;
		}
	}*/
	
	
	file.open("flow_statistics.txt");
	file << "u0: " << region.u0 << std::endl;
	file << "turbulence intensity: " << TI1 << std::endl;
	file.close();


	for (size_t i{ 0 }; i < iters; i++) {
		region.increment_eddies();
		region.print_flucts(i, "MRSEM_test");
		std::cout << i << std::endl;
	}


	
	/*
	for (size_t i{ 0 }; i < region.eddies.size; i++) {
		std::cout << std::setw(6) << region.eddies(i).position(0) << ", ";
		std::cout << std::setw(6) << region.eddies(i).position(1) << ", ";
		std::cout << std::setw(6) << region.eddies(i).position(2) << "| ";
		std::cout << std::setw(6) << region.eddies(i).radius << std::endl;
	}
	*/
	//region.delta = 0.004;



}