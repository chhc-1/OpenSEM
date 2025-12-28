#pragma once

#include "Array.h"
#include "MRSEM_eddy.h"
#include "MRSEM_subregion.h"
#include "MRSEM_region.h"

#define e 2.7182818284590452353

double gaussian(const double& x) {
	//return exp(-0.5 * x*x) * 1 / 0.8641898707968644; // factor comes from evaluation of gaussian squared from -1 to 1 - see jupyter notebook
	return exp(-0.5 * x * x);
}
double mgaussian(const double& x) {
	return -exp(-0.5 * x * x);// *1 / 0.8641898707968644;
}

// x gaus may have issues
double xgaussian(const double& x) {
	return 1 / 0.4352842126938356 * std::abs(x) * exp(-0.5 * x*x);
}

// H function most likely has issues - all regions with H have issues currently
double H(const double& x){
	//return 1 / 0.7071067811865476 * (1 - cos(2 * M_PI * std::abs(x)) / sqrt(0.214));// (2 * M_PI * sqrt(0.214));
	return cos(M_PI / 2 * x);// *sqrt(1.6);// *sqrt(2);
}

double zerofnptr(const double& x) {
	return 0;
}

double binaryfnptr(const double& x) {
	return 1;
}


// need to scale the functions separately for each region, since eg a larger eddy size in the z direction gives higher probability of non-zero value -> increases value
// Increase in value is likely reflected by spike in produced RST  
// May not be necessary since number of eddies is also scaled by this
// evaluate the scaling factor numerically, also need to consider the effect of different time lengths (lt)
// also consider taking the product of the scaling functions and producing the scaling factor using that, although this should not be necessary if number eddies
// balances the eddy size/radius appropriately

// currently output RST is independent of input eddy radius (assuming same ratio of eddy radii), which is as expected
// values however are of order 10^2 off...

MRSEM_region TBL_MRSEM(const size_t& ny, const size_t& nz, const double& ymax, const double& zmax, const double& _u0, const double& _dt,
	const double& delta, const double& utau, const double& nu) {


	shape_fn gauss1(gaussian);
	shape_fn mgauss1(mgaussian);
	shape_fn xgauss1(xgaussian);
	shape_fn H1(H);
	shape_fn zerofn(zerofnptr);
	shape_fn binaryfn(binaryfnptr);

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

	size_t nregions = 3;//4;//5; 
	// currently using: 1 near wall region, 1 hairpin region (tails), 1 outer region

	MRSEM_region _MRSEM = MRSEM_region(y, z, _u0, nregions);

	Array<double> y_temp;
	Array<double> z_temp;

	y_temp.resize({ny, nz});
	z_temp.resize({ ny, nz });
		
	y_temp.set(0);
	z_temp.set(0);

	std::cout << "length scale: "  << l << std::endl;

	// near wall region
	{
		double lx1 = 100 * l; // 60 * l; //
		double ly1 = 20 * l; // 60 * l; //
		double lz1 = 60 * l; // 60 * l; //
		double c1 = 15 * utau;

		//std::cout << _MRSEM.y_inlet.shape[0] << ", " << _MRSEM.z_inlet.shape[1] << std::endl;
		_MRSEM.regions(0) = MRSEM_subregion(c1, _dt, xtemp, _MRSEM.y_inlet, _MRSEM.z_inlet, rad_temp, delta, 
			lx1, ly1, lz1, 20 * l, 60 * l);
		//std::cout << _MRSEM.y_inlet.shape[0] << ", " << _MRSEM.z_inlet.shape[1] << std::endl;
		_MRSEM.regions(0).u0 = c1;
		_MRSEM.regions(0).dt = _dt;
		_MRSEM.regions(0).delta = delta;
		_MRSEM.regions(0).lx = lx1;
		//_MRSEM.regions(0).ly= /

		_MRSEM.regions(0).shape_t[0] = gauss1;//binaryfn; //gauss1;
		_MRSEM.regions(0).shape_t[1] = gauss1;//binaryfn; //gauss1;
		_MRSEM.regions(0).shape_t[2] = H1;//binaryfn; //H1;//gauss1;//H1;

		_MRSEM.regions(0).shape_y[0] = gauss1;//binaryfn; //gauss1;
		_MRSEM.regions(0).shape_y[1] = gauss1;//binaryfn; //mgauss1;
		_MRSEM.regions(0).shape_y[2] = H1;//binaryfn; //H1;//gauss1;//H1;

		_MRSEM.regions(0).shape_z[0] = gauss1;//binaryfn; //gauss1;
		_MRSEM.regions(0).shape_z[1] = gauss1;//binaryfn; //;
		_MRSEM.regions(0).shape_z[2] = H1;//binaryfn; //H1;//gauss1;//H1;

		double temp = _MRSEM.regions(0).shape_t[0](1, 0.5, 1);
		//std::cout << "shape eval: " << _MRSEM.regions(0).shape_t[0](1, 0.5, 1) << ", " << gauss1(1, 0.5, 1) << std::endl;

		_MRSEM.regions(0).vf_scaling_factor *= (2* 0.8641898707968644 )*sqrt(2.2135805);

		// 20 < y+ < 60
		_MRSEM.regions(0).y_rand = std::uniform_real_distribution<double>(20 * l, 60 * l); 

		//std::cout << _MRSEM.regions(0).d_max << ", " << _MRSEM.y_inlet.shape[0] << ", " << _MRSEM.z_inlet.shape[1] << std::endl;
	}

	std::cout << "created near wall region. eddies:" << _MRSEM.regions(0).eddies.size << std::endl;
	// hair pin vortices / log layer
	// hairpin legs
	{
		double lx2 = 120 * l; //60 * l; //
		double ly2 = 60 * l; //60 * l; //30 * l; //
		double lz2 = 60 * l; //60 * l; // 30 * l; //
		double c2 = 15 * utau;
		//std::cout << _MRSEM.y_inlet.shape[0] << ", " << _MRSEM.z_inlet.shape[1] << std::endl;
		//_MRSEM.regions(1) = MRSEM_subregion(c2, _dt, xtemp, _MRSEM.y_inlet, _MRSEM.z_inlet, rad_temp, delta, lx2, ly2, lz2,
		//	60 * l, 0.4 * delta);
		_MRSEM.regions(1) = MRSEM_subregion(c2, _dt, xtemp, _MRSEM.y_inlet, _MRSEM.z_inlet, rad_temp, delta, lx2, ly2, lz2,
				60 * l, 0.5 * delta);

		_MRSEM.regions(1).shape_t[0] = gauss1; //binaryfn; //gauss1;
		_MRSEM.regions(1).shape_t[1] = gauss1; //binaryfn; //gauss1;
		_MRSEM.regions(1).shape_t[2] = H1;//binaryfn; //H1;//gauss1;//H1;

		_MRSEM.regions(1).shape_y[0] = gauss1; //binaryfn; //gauss1;
		_MRSEM.regions(1).shape_y[1] = gauss1; //binaryfn; //mgauss1;
		_MRSEM.regions(1).shape_y[2] = H1;//binaryfn; //H1;//gauss1;//H1;

		_MRSEM.regions(1).shape_z[0] = gauss1; //binaryfn; //gauss1;
		_MRSEM.regions(1).shape_z[1] = gauss1; //binaryfn; //gauss1;
		_MRSEM.regions(1).shape_z[2] = H1;//binaryfn; //H1;//gauss1;//H1;
		
		_MRSEM.regions(1).vf_scaling_factor *= (2* 0.8641898707968644)*sqrt(2.2135805);

		// 60 < y+ < 0.4deltaplus
		//_MRSEM.regions(1).y_rand = std::uniform_real_distribution<double>(60 * l, 0.4 * delta);
		_MRSEM.regions(1).y_rand = std::uniform_real_distribution<double>(60 * l, 0.5 * delta);
	}

	std::cout << "created hairpin tail region. eddies:" << _MRSEM.regions(1).eddies.size << std::endl;
	/*
	// hairpin heads
	{
		double lx3 = 60 * l; //60 * l;
		double ly3 = 60 * l; //60 * l;
		double lz3 = 120 * l; //120 * l; //120 * l;//80 * l;//120 * l;
		double c3 = 15 * utau;

		_MRSEM.regions(2) = MRSEM_subregion(c3, _dt, xtemp, _MRSEM.y_inlet, _MRSEM.z_inlet, rad_temp, delta, lx3, ly3, lz3,
			0.4 * delta, 0.5 * delta);

		_MRSEM.regions(2).shape_t[0] = zerofn; //mgauss1;
		_MRSEM.regions(2).shape_t[1] = zerofn; //H1;//gauss1;//H1;
		_MRSEM.regions(2).shape_t[2] = zerofn; //gauss1;

		_MRSEM.regions(2).shape_y[0] = zerofn; //xgauss1;//gauss1;//xgauss1;
		_MRSEM.regions(2).shape_y[1] = zerofn; //gauss1;
		_MRSEM.regions(2).shape_y[2] = zerofn; //gauss1;

		_MRSEM.regions(2).shape_z[0] = zerofn; //gauss1;
		_MRSEM.regions(2).shape_z[1] = zerofn; //gauss1;
		_MRSEM.regions(2).shape_z[2] = zerofn; //xgauss1;//gauss1;//xgauss1;

		// 0.4deltaplus < y+ < 0.5deltaplus
		_MRSEM.regions(2).y_rand = std::uniform_real_distribution<double>(0.4 * delta, 0.5 * delta);
	}*/

	//std::cout << "created hairpin head region. eddies:" << _MRSEM.regions(2).eddies.size << std::endl;

	// outer layers
	// outer region A
	{
		double lx4 = 0.1 * delta;
		double ly4 = 0.1 * delta;
		double lz4 = 0.1 * delta;
		double c4 = 0.8 * _u0;

		//_MRSEM.regions(3) = MRSEM_subregion(c4, _dt, xtemp, _MRSEM.y_inlet, _MRSEM.z_inlet, rad_temp, delta, lx4, ly4, lz4,
		//	0.5 * delta, 0.8 * delta);
		_MRSEM.regions(2) = MRSEM_subregion(c4, _dt, xtemp, _MRSEM.y_inlet, _MRSEM.z_inlet, rad_temp, delta, lx4, ly4, lz4,
			0.5 * delta, 1.0 * delta);

		_MRSEM.regions(2).shape_t[0] = gauss1;//binaryfn; //gauss1;
		_MRSEM.regions(2).shape_t[1] = gauss1;//binaryfn; //gauss1;
		_MRSEM.regions(2).shape_t[2] = gauss1;//binaryfn; //gauss1;

		_MRSEM.regions(2).shape_y[0] = gauss1;//binaryfn; //gauss1;
		_MRSEM.regions(2).shape_y[1] = gauss1;//binaryfn; //gauss1;
		_MRSEM.regions(2).shape_y[2] = gauss1;//binaryfn; //gauss1;

		_MRSEM.regions(2).shape_z[0] = gauss1;//binaryfn; //gauss1;
		_MRSEM.regions(2).shape_z[1] = gauss1;//binaryfn; //gauss1;
		_MRSEM.regions(2).shape_z[2] = gauss1;//binaryfn; //gauss1;

		_MRSEM.regions(2).vf_scaling_factor *= (2*0.8641898707968644 )*sqrt(1.680062502);

		//_MRSEM.regions(3).y_rand = std::uniform_real_distribution<double>(0.5 * delta, 0.8 * delta);
		_MRSEM.regions(2).y_rand = std::uniform_real_distribution<double>(0.5 * delta, 1.0 * delta);
	}
	std::cout << "created outer region A. eddies:" << _MRSEM.regions(2).eddies.size << std::endl;
	/*
	// outer region B
	{
		double lx5 = 0.15 * delta;
		double ly5 = 0.15 * delta;
		double lz5 = 0.15 * delta;
		double c5 = 0.8 * _u0;

		_MRSEM.regions(4) = MRSEM_subregion(c5, _dt, xtemp, _MRSEM.y_inlet, _MRSEM.z_inlet, rad_temp, delta, lx5, ly5, lz5,
			0.8 * delta, y.max());

		_MRSEM.regions(4).shape_t[0] = binaryfn; //gauss1;
		_MRSEM.regions(4).shape_t[1] = binaryfn; //gauss1;
		_MRSEM.regions(4).shape_t[2] = binaryfn; //gauss1;

		_MRSEM.regions(4).shape_y[0] = binaryfn; //gauss1;
		_MRSEM.regions(4).shape_y[1] = binaryfn; //gauss1;
		_MRSEM.regions(4).shape_y[2] = binaryfn; //gauss1;

		_MRSEM.regions(4).shape_z[0] = binaryfn; //gauss1;
		_MRSEM.regions(4).shape_z[1] = binaryfn; //gauss1;
		_MRSEM.regions(4).shape_z[2] = binaryfn; //gauss1;

		_MRSEM.regions(4).y_rand = std::uniform_real_distribution<double>(0.8 * delta, y.max());
	}
	std::cout << "created outer region B. eddies:" << _MRSEM.regions(4).eddies.size << std::endl;
	*/
	return _MRSEM;
}