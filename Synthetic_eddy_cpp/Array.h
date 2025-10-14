#pragma once

//#include <Vector.h>
#include <malloc.h>
#include <initializer_list>
#include <cassert>
#include <cmath>
#include <exception>
#include <iostream>

template <typename T>
struct Array { // assume array is 2D, vectors have 1 length in 2nd dimension
public:
    size_t size{ 0 };
    size_t ndims{ 0 };
    size_t* shape = nullptr; // array 
    T* array1 = nullptr; // array

    Array() {

    }
    /*
    Array<T> operator=(const std::initializer_list<T> init_values) {

    }
    */

    Array(std::initializer_list<size_t> list) {
        shape = new size_t[list.size()];
        ndims = list.size();
        size = 1;
        for (size_t i{ 0 }; i < list.size(); i++) {
            shape[i] = *(list.begin() + i);
            size *= shape[i];
        }
        array1 = new T[size];
    }
    
    /*
    ~Array(){
        free(shape);
        free(array1);
    }*/

    /*
    Array(const size_t _shape[]) {
        shape = new size_t[sizeof(_shape) / sizeof(T)];
        size = 1;
        for (size_t i{ 0 }; i < sizeof(_shape) / sizeof(T); i++) {
            shape[i] = _shape[i];
            size *= _shape[i];
        }
        //ndims = sizeof(_shape) / sizeof(T);

        array1 = (T*)malloc(size * sizeof(T));
    }

    Array(const size_t* _shape[]) {
        shape = new size_t[sizeof(*_shape) / sizeof(T)];
        //ndims = sizeof(_shape) / sizeof(T);
        for (size_t i{ 0 }; i < sizeof(*_shape) / sizeof(T); i++) {
            shape[i] = _shape[i];
            size *= _shape[i];
        }
        array1 = (T*)malloc(size * sizeof(T));
    }
    */

    /*
    void resize(const size_t* _shape[]) {
        delete shape;
        shape = new size_t[sizeof(*_shape) / sizeof(size_t)];
        //ndims = sizeof(*_shape) / sizeof(T);
        size = 1;
        for (size_t i{ 0 }; i < sizeof(*_shape) / sizeof(size_t); i++) {
            shape[i] = _shape[i];
            size *= shape[i];
        }
        delete(array1);
        array1 = (T*)malloc(size * sizeof(T));
    }
    */
    void resize(const size_t _shape[], const size_t& _ndims) {
        //std::cout << sizeof(*_shape) / sizeof(size_t) << std::endl;
        delete shape;
        //shape = new size_t[sizeof(*_shape) / sizeof(size_t)];
        shape = new size_t[_ndims];
        ndims = _ndims;
        //ndims = sizeof(*_shape) / sizeof(T);
        size = 1;
        for (size_t i{ 0 }; i < _ndims; i++) {
            shape[i] = _shape[i];
            size *= shape[i];
        }
        delete(array1);
        array1 = (T*)malloc(size * sizeof(T));
    }

    void resize(const std::initializer_list<size_t> list) {

        delete shape;
        shape = new size_t[list.size()];
        ndims = list.size();
        size = 1;
        for (size_t i{ 0 }; i < ndims; i++) {
            shape[i] = *(list.begin() + i);
            size *= shape[i];
        }
        delete(array1);
        array1 = (T*)malloc(size * sizeof(T));
    }

    void set(const T& inp_value) {
        for (size_t i{ 0 }; i < size; i++) {
            array1[i] = inp_value;
        }
    }

    // add block indexing?

    T& get(const size_t& idx) {
        assert(idx < size);
        return array1[idx];
    }

    T& get(const size_t& idx) const {
        assert(idx < size);
        return array1[idx];
    }

    T& get(const size_t& s1, const size_t& s2) {
        //assert(ndims > 1);
        assert(s1 < shape[0]);
        assert(s2 < shape[1]);
        //assert(shape[1] * s1 + s2 < size);
        return array1[shape[1] * s1 + s2];
    }

    const T& get(const size_t& s1, const size_t& s2) const {
        //assert(ndims > 1);
        assert(shape[1] * s1 + s2 < size);
        return array1[shape[1] * s1 + s2];
    }

    T& get(size_t& s1, size_t& s2) {
        //assert(ndims > 1);
        assert(shape[1] * s1 + s2 < size);
        return array1[shape[1] * s1 + s2];
    }

    T& operator()(const size_t& idx) {
        return get(idx);
    }

    T& operator()(const size_t& idx) const {
        return get(idx);
    }

    T& operator()(const size_t& s1, const size_t& s2) {
        return get(s1, s2);
    }

    const T& operator()(const size_t& s1, const size_t& s2) const {
        return get(s1, s2);
    }

    T& operator()(size_t& s1, size_t& s2) {
        return get(s1, s2);
    }

    const T min() { // look for min in entire array
        size_t idx_1D{ 0 };
        T value = array1[0];
        for (size_t i{ 1 }; i < size; i++) {
            if (array1[i] < value) {
                value = array1[i];
                idx_1D = i;
            }
        }
        size_t* idx = new size_t[2];
        auto dv = std::div(idx_1D, (int)shape[1]);
        idx[0] = dv.quot;
        idx[1] = dv.rem;
        return array1[idx[0] * shape[1] + idx[1]];
    }


    const T min(const size_t& col) { // look for min in given row
        size_t idx{ 0 };
        T value = array1[col];
        for (size_t i{ 1 }; i < shape[0]; i++) {
            if (array1[col + i * shape[1]] < value) {
                value = array1[col + i * shape[1]];
                idx = i;
            }
        }
        return array1[idx];
    }

    const T max() { // look for max in entire array
        size_t idx_1D{ 0 };
        T value = array1[0];
        for (size_t i{ 1 }; i < size; i++) {
            if (array1[i] > value) {
                value = array1[i];
                idx_1D = i;
            }
        }
        size_t* idx = new size_t[2];
        auto dv = std::div(idx_1D, (int)shape[1]);
        idx[0] = dv.quot;
        idx[1] = dv.rem;
        return array1[idx[0] * shape[1] + idx[1]];
    }

    const T max(const size_t& col) { // look for min in given row
        size_t idx{ 0 };
        T value = array1[col];
        for (size_t i{ 1 }; i < shape[0]; i++) {
            if (array1[col + i * shape[1]] > value) {
                value = array1[col + i * shape[1]];
                idx = i;
            }
        }
        return array1[idx];
    }

    size_t* where(const T& value) { // look for matching value in entire array
        //size_t idx_1D{ 0 };
        for (size_t i{ 0 }; i < size; i++) {
            if (value == array1[i]) {
                auto dv = std::div(int(i), (int)shape[1]);
                size_t* idx = new size_t[2];
                idx[0] = dv.quot;
                idx[1] = dv.rem;
                return idx;
                break;
                //idx_1D = i;
                //break;
            }
        }
        throw std::exception("Value not found");
    }

    const size_t where(const T& value, const size_t& col) { // look for matching value in given col
        for (size_t i{ 0 }; i < shape[0]; i++) {
            if (value == array1[col + i * shape[1]]) {
                return i;
                break;
            }
        }
        throw std::exception("Value not found");
        //throw std::nested_exception("Did not find value");
    }

    void pop_row(const size_t& row_num) {
        for (size_t i{ row_num }; i < shape[0] - 1; i++) {
            for (size_t j{ 0 }; j < shape[1]; j++) {
                array1[i * shape[1] + j] = array1[(i + 1) * shape[1] + j];
            }
        }
        for (size_t j{ 0 }; j < shape[1]; j++) {
            array1[(shape[0] - 1) * shape[1] + j] = 0;
        }
    }

    void pop_row(const size_t& row_num, size_t& cutoff) {
        for (size_t i{ row_num }; i < cutoff - 1; i++) {
            for (size_t j{ 0 }; j < shape[1]; j++) {
                array1[i * shape[1] + j] = array1[(i + 1) * shape[1] + j];
            }
        }
        for (size_t j{ 0 }; j < shape[1]; j++) {
            array1[cutoff * shape[1] + j] = 0;
        }
        cutoff -= 1;
    }

    void pop_col(const size_t& col_num, size_t& cutoff) {
        for (size_t j{ col_num }; j < cutoff - 1; j++) {
            for (size_t i{ 0 }; i < shape[1]; i++) {
                array1[i * shape[0] + j] = array1[i * shape[0] + j + 1];
            }
        }
        for (size_t i{ 0 }; i < shape[0]; i++) {
            array1[i * shape[0] + cutoff] = 0;
        }
        cutoff -= 1;
    }

    const bool& is_in(const T& value) const {
        for (size_t i{ 0 }; i < size; i++) {
            if (array1[i] == value) {
                return true;
            }
        }
        return false;
    }

    const bool& in_col(const T& value, const size_t& col_num) const {
        for (size_t i{ 0 }; i < shape[0]; i++) {
            if (array1[i * shape[1] + col_num] == value) {
                return true;
            }
        }
        return false;
    }

    const bool& in_row(const T& value, const size_t& row_num) const {
        for (size_t i{ 0 }; i < shape[1]; i++) {
            if (array1[row_num * shape[1] + i] == value) {
                return true;
            }
        }
        return false;
    }

    void LinRange(const T& start, const T& end, const size_t& n_points) {
        double increment = (end - start) / (n_points - 1);
        resize({ n_points });
        array1[0] = start;
        for (size_t i{ 1 }; i < n_points; i++) {
            array1[i] = start + increment * i;
        }
    }

    void mul(const T& mult) {
        for (size_t a{ 0 }; a < size; a++) {
            array1[a] *= mult;
        }
    }

    /*

    void LinRange(const T& start, const T& end) {
        double increment = (end - start) / (n_points - 1);
        for (size_t i{ 0 }; i < n_points; i++) {

        }
    }

    void LinRange(const T& start, const T& end, const size_t& n_points) {
        double increment = (end - start) / (n_points - 1);
        for (size_t i{ 0 }; i < n_points; i++) {

        }
    }
    */
    /*
    void mul() { // in place multiplication over array

    }
    */
    /*private:



    public:
    */

};

