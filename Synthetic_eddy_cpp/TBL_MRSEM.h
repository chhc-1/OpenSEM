#pragma once

#include "Array.h"
#include "MRSEM_eddy.h"
#include "MRSEM_subregion.h"
#include "MRSEM_region.h"

double gaussian(const double& x) {
	return 1 / sqrt(2 * M_PI) * exp(-0.5 * x);
}
double mgaussian(const double& x) {
	return -1 / sqrt(2 * M_PI) * exp(-0.5 * x);
}
double xgaussian(const double& x) {
	return  x / sqrt(2 * M_PI) * exp(-0.5 * x);
}
double H(const double& x){
	return 1 - cos(2 * M_PI * x) / (2 * M_PI * x * sqrt(0.214));
}


MRSEM_region TBL_MRSEM(const size_t& ny, const size_t& nz, const double& ymax, const double& zmax, const double& _u0, const double& _dt,
	const double& delta, const double& utau, const double& nu) {


	shape_fn gauss1(gaussian);
	shape_fn mgauss1(mgaussian);
	shape_fn xgauss1(xgaussian);
	shape_fn H1(H);

	//std::cout << gauss1(0.5) << ", " << gauss1(1, 0.5, 1) << std::endl;
	

	double l = nu / utau;
	double lx;
	double ly;
	double lz;
	double deltaplus = delta / l;

	double xtemp = 1;
	double rad_temp = 0;
	Array<double> y;
	Array<double> yplus;

	size_t ystart[5];
	size_t yend[5];

	Array<double> z;
	
	std::cout << ymax << ", " << zmax << std::endl;

	y.LinRange(0, ymax, ny);
	z.LinRange(0, zmax, nz);
	yplus.resize({ y.size });
	for (size_t i{ 0 }; i < y.size; i++) {
		yplus(i) = y(i) / l;
	}

	size_t nregions = 5;

	MRSEM_region _MRSEM = MRSEM_region(y, z, _u0, nregions);

	Array<double> y_temp;
	Array<double> z_temp;

	y_temp.resize({ny, nz});
	z_temp.resize({ ny, nz });
		
	y_temp.set(0);
	z_temp.set(0);

	// near wall region
	{
		double lx1 = 100 * l;
		double ly1 = 20 * l;
		double lz1 = 60 * l;
		double c1 = 15 * utau;

		//std::cout << _MRSEM.y_inlet.shape[0] << ", " << _MRSEM.z_inlet.shape[1] << std::endl;
		_MRSEM.regions(0) = MRSEM_subregion(c1, _dt, xtemp, _MRSEM.y_inlet, _MRSEM.z_inlet, rad_temp, delta, lx1, ly1, lz1);
		//std::cout << _MRSEM.y_inlet.shape[0] << ", " << _MRSEM.z_inlet.shape[1] << std::endl;
		_MRSEM.regions(0).u0 = c1;
		_MRSEM.regions(0).dt = _dt;
		_MRSEM.regions(0).delta = delta;
		_MRSEM.regions(0).lx = lx1;
		//_MRSEM.regions(0).ly= /

		_MRSEM.regions(0).shape_t[0] = gauss1;
		_MRSEM.regions(0).shape_t[1] = gauss1;
		_MRSEM.regions(0).shape_t[2] = H1;

		_MRSEM.regions(0).shape_y[0] = gauss1;
		_MRSEM.regions(0).shape_y[1] = mgauss1;
		_MRSEM.regions(0).shape_y[2] = H1;

		_MRSEM.regions(0).shape_z[0] = gauss1;
		_MRSEM.regions(0).shape_z[1] = gauss1;
		_MRSEM.regions(0).shape_z[2] = H1;

		double temp = _MRSEM.regions(0).shape_t[0](1, 0.5, 1);
		//std::cout << "shape eval: " << _MRSEM.regions(0).shape_t[0](1, 0.5, 1) << ", " << gauss1(1, 0.5, 1) << std::endl;

		// 20 < y+ < 60
		_MRSEM.regions(0).y_rand = std::uniform_real_distribution<double>(20 * l, 60 * l); 

		//std::cout << _MRSEM.regions(0).d_max << ", " << _MRSEM.y_inlet.shape[0] << ", " << _MRSEM.z_inlet.shape[1] << std::endl;
	}

	std::cout << "created near wall region. eddies:" << _MRSEM.regions(0).eddies.size << std::endl;
	// hair pin vortices / log layer
	// hairpin legs
	{
		double lx2 = 120 * l;
		double ly2 = 60 * l;
		double lz2 = 60 * l;
		double c2 = 15 * utau;
		//std::cout << _MRSEM.y_inlet.shape[0] << ", " << _MRSEM.z_inlet.shape[1] << std::endl;
		_MRSEM.regions(1) = MRSEM_subregion(c2, _dt, xtemp, _MRSEM.y_inlet, _MRSEM.z_inlet, rad_temp, delta, lx2, ly2, lz2);

		_MRSEM.regions(1).shape_t[0] = gauss1;
		_MRSEM.regions(1).shape_t[1] = gauss1;
		_MRSEM.regions(1).shape_t[2] = H1;

		_MRSEM.regions(1).shape_y[0] = gauss1;
		_MRSEM.regions(1).shape_y[1] = mgauss1;
		_MRSEM.regions(1).shape_y[2] = H1;

		_MRSEM.regions(1).shape_z[0] = gauss1;
		_MRSEM.regions(1).shape_z[1] = gauss1;
		_MRSEM.regions(1).shape_z[2] = H1;
		
		// 60 < y+ < 0.4deltaplus
		_MRSEM.regions(1).y_rand = std::uniform_real_distribution<double>(60 * l, 0.4 * delta);
	}

	std::cout << "created hairpin tail region. eddies:" << _MRSEM.regions(1).eddies.size << std::endl;

	// hairpin heads
	{
		double lx3 = 60 * l;
		double ly3 = 60 * l;
		double lz3 = 120 * l;
		double c3 = 15 * utau;

		_MRSEM.regions(2) = MRSEM_subregion(c3, _dt, xtemp, _MRSEM.y_inlet, _MRSEM.z_inlet, rad_temp, delta, lx3, ly3, lz3);

		_MRSEM.regions(2).shape_t[0] = mgauss1;
		_MRSEM.regions(2).shape_t[1] = H1;
		_MRSEM.regions(2).shape_t[2] = gauss1;

		_MRSEM.regions(2).shape_y[0] = xgauss1;
		_MRSEM.regions(2).shape_y[1] = gauss1;
		_MRSEM.regions(2).shape_y[2] = gauss1;

		_MRSEM.regions(2).shape_z[0] = gauss1;
		_MRSEM.regions(2).shape_z[1] = gauss1;
		_MRSEM.regions(2).shape_z[2] = xgauss1;

		// 0.4deltaplus < y+ < 0.5deltaplus
		_MRSEM.regions(2).y_rand = std::uniform_real_distribution<double>(0.4 * delta, 0.5 * delta);
	}

	std::cout << "created hairpin head region. eddies:" << _MRSEM.regions(2).eddies.size << std::endl;

	// outer layers
	// outer region A
	{
		double lx4 = 0.1 * delta;
		double ly4 = 0.1 * delta;
		double lz4 = 0.1 * delta;
		double c4 = 0.8 * _u0;

		_MRSEM.regions(3) = MRSEM_subregion(c4, _dt, xtemp, _MRSEM.y_inlet, _MRSEM.z_inlet, rad_temp, delta, lx4, ly4, lz4);

		_MRSEM.regions(3).shape_t[0] = gauss1;
		_MRSEM.regions(3).shape_t[1] = gauss1;
		_MRSEM.regions(3).shape_t[2] = gauss1;

		_MRSEM.regions(3).shape_y[0] = gauss1;
		_MRSEM.regions(3).shape_y[1] = gauss1;
		_MRSEM.regions(3).shape_y[2] = gauss1;

		_MRSEM.regions(3).shape_z[0] = gauss1;
		_MRSEM.regions(3).shape_z[1] = gauss1;
		_MRSEM.regions(3).shape_z[2] = gauss1;

		_MRSEM.regions(3).y_rand = std::uniform_real_distribution<double>(0.5 * delta, 0.8 * delta);
	}
	std::cout << "created outer region A. eddies:" << _MRSEM.regions(3).eddies.size << std::endl;

	// outer region B
	{
		double lx5 = 0.15 * delta;
		double ly5 = 0.15 * delta;
		double lz5 = 0.15 * delta;
		double c5 = 0.8 * _u0;

		_MRSEM.regions(4) = MRSEM_subregion(c5, _dt, xtemp, _MRSEM.y_inlet, _MRSEM.z_inlet, rad_temp, delta, lx5, ly5, lz5);

		_MRSEM.regions(4).shape_t[0] = gauss1;
		_MRSEM.regions(4).shape_t[1] = gauss1;
		_MRSEM.regions(4).shape_t[2] = gauss1;

		_MRSEM.regions(4).shape_y[0] = gauss1;
		_MRSEM.regions(4).shape_y[1] = gauss1;
		_MRSEM.regions(4).shape_y[2] = gauss1;

		_MRSEM.regions(4).shape_z[0] = gauss1;
		_MRSEM.regions(4).shape_z[1] = gauss1;
		_MRSEM.regions(4).shape_z[2] = gauss1;

		_MRSEM.regions(4).y_rand = std::uniform_real_distribution<double>(0.8 * delta, y.max());
	}
	std::cout << "created outer region B. eddies:" << _MRSEM.regions(4).eddies.size << std::endl;

	return _MRSEM;
}