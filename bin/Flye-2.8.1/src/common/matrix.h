//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <cstdlib>
#include <stdexcept>
#include <cassert>
#include <vector>

template <class T>
class Matrix
{
public:
	Matrix(): _rows(0), _cols(0), _data(nullptr) {}
	Matrix(const Matrix& other):
		Matrix(other._rows, other._cols)
	{
		for (size_t i = 0; i < _rows; ++i)
			for (size_t j = 0; j < _cols; ++j)
				this->at(i, j) = other.at(i, j);
	}
	Matrix(Matrix&& other):
		Matrix()
	{
		std::swap(_cols, other._cols);
		std::swap(_rows, other._rows);
		std::swap(_data, other._data);
	}
	Matrix& operator=(Matrix && other)
	{
		std::swap(_cols, other._cols);
		std::swap(_rows, other._rows);
		std::swap(_data, other._data);
		return *this;
	}
	Matrix& operator=(const Matrix& other)
	{
		Matrix temp(other);
		std::swap(_cols, temp._cols);
		std::swap(_rows, temp._rows);
		std::swap(_data, temp._data);
		return *this;
	}

	Matrix(size_t rows, size_t cols, T val = 0):
		_rows(rows), _cols(cols)
	{
		if (!rows || !cols)
			throw std::runtime_error("Zero matrix dimension");
		_data = new T[rows * cols];
		for (size_t i = 0; i < rows * cols; ++i) _data[i] = val;
	}
	~Matrix()
	{
		if (_data) delete[] _data;
	}

	T& at(size_t row, size_t col) {return _data[row * _cols + col];}
	const T& at(size_t row, size_t col) const {return _data[row * _cols + col];}
	size_t nrows() const {return _rows;}
	size_t ncols() const {return _cols;}

private:
	size_t _rows;
	size_t _cols;
	T* _data;
};
