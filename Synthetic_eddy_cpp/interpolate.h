#pragma once

// MUST FIX

#include "Array.h"
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <iostream>

template <typename T>
struct interpolator {
	Array<T> x;
	Array<T> y;

	interpolator(){

	}
	interpolator(Array<T> _x, Array<T> _y) {
		x = _x;
		y = _y;
	}

	const T interpolate(const T& value) {
		size_t min_idx{ 0 };
		T err;
		T err_prev;

		assert(value > x(0));
		assert(value < x(x.size - 1));

		err_prev = value - x(0);

		for (size_t i{ 1 }; i < x.size; i++) {
			err = value - x(i);
			if (std::abs(err) < std::abs(err_prev) && err > 0) {
				err_prev = err;
				min_idx = i;
			}
			if (err < 0) {
				//std::cout << min_idx << std::endl;
				break;
			}
		}
		return (y(min_idx+1) - y(min_idx)) / (x(min_idx+1) - x(min_idx)) * (value-x(min_idx)) + y(min_idx);
	}

	const T operator()(const T& value) {
		return interpolate(value);
	}

};