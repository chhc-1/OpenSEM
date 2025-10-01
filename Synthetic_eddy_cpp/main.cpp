#include "Array.h"
#include "Eddy.h"
#include "oSEM_eddy.h"
#include "oSEM_region.h"
#include "ISEM1_eddy.h"
#include "ISEM1_region.h"
#include "DFSEM_eddy.h"
#include "DFSEM_region.h"

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
};



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

	double delta = 0.004;

	double x_pos = 0.0;
	double y_min = 0.0;
	double y_max = 1.4 * delta;
	double z_min = 0.0;
	double z_max = 0.03;

	size_t n_y = 50; // 200;
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

	/*
	Array<double> temp1;
	temp1.resize(y_plane.shape, y_plane.ndims);
	std::cout << sizeof(y_plane.shape[0]) << std::endl;
	std::cout << sizeof(&(y_plane.shape)) / sizeof(size_t) << std::endl;
	std::cout << temp1.shape[0] << ", " << temp1.shape[1] << std::endl;
	*/
	//size_t arr[] = {1 ,3 ,5};
	//size_t (*ptr)[3] = &(arr);

	double u0 = 150; // 100;
	double dt = 0.000002;

	ISEM1_region region = ISEM1_region(u0, dt, x_pos, y_plane, z_plane, rep_rad1, max_rad1, delta);
	//oSEM_region region = oSEM_region(u0, dt, x_pos, y_plane, z_plane, max_rad1, delta);
	//DFSEM_region region = DFSEM_region(u0, dt, x_pos, y_plane, z_plane, min_rad1, rep_rad1, max_rad1, delta);

	double TI1 = 0.01;

	region.set_HIT_RST(TI1);

	//_tracker.print_allocs();
	//region.increment_eddies();
	//_tracker.print_allocs();
	//region.increment_eddy(region.eddies(0));
	//_tracker.print_allocs();

	//std::uniform_int_distribution<int> dist3 = std::uniform_int_distribution<int>(0, 1);
	//std::mt19937 urng = std::mt19937(1729);

	//int test1;

	
	//test1 = dist3(urng);
	//_tracker.print_allocs();

	std::cout << region.eddies.size << std::endl;

	//std::cout << region.a11(0) << std::endl;

	//std::cout << region.vf_scaling_factor << std::endl;
	//std::cout << region.eddies(0).shape_scaling_factor << std::endl;
	//std::cout << region.vol << ", " << pow(region.vol, 0.5) << ", " << region.vol_sqrt << std::endl;
	//std::cout << region.vol / pow(rep_rad1, 3) << std::endl;
	//std::cout << region.vf_scaling_factor * region.eddies(0).shape_scaling_factor << std::endl;

	size_t iters{ 1000 };

	std::ofstream file;

	file.open("flow_statistics.txt");
	file << "u0: " << region.u0 << std::endl;
	file << "turbulence intensity: " << TI1 << std::endl;
	file.close();

	for (size_t i{ 0 }; i < iters; i++) {
		region.increment_eddies();
		region.print_flucts(i, "ISEM_test");
		std::cout << i << std::endl;
	}


	/*
	for (size_t i{ 0 }; i < region.a11.size; i++) {
		std::cout << region.a11(i) << ", ";
		std::cout << region.a21(i) << ", ";
		std::cout << region.a22(i) << ", ";
		std::cout << region.a31(i) << ", ";
		std::cout << region.a32(i) << ", ";
		std::cout << region.a33(i) << std::endl;
	}
	*/
	
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