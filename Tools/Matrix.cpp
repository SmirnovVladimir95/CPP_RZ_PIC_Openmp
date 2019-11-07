#include "Matrix.h"

Matrix::Matrix(size_t rows, size_t columns) : _rows(rows), _columns(columns) {
    data.resize(rows*columns);
}

size_t Matrix::rows() const { return _rows; }

size_t Matrix::columns() const { return _columns; }

type_double* Matrix::data_ptr() { return data.data(); }

type_double& Matrix::operator()(size_t row, size_t column) {
    if (row >= _rows or column >= _columns) {
        cout << "Error, index (" << row << " " << column << ") out of range (return 1st element)" << endl;
        return data[0];
    }
    return data[row * _columns + column];
}

void Matrix::print() {
    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _columns; j++) {
            cout << this->operator()(i, j) << " ";
        }
        cout << endl;
    }
}

void Matrix::fill(const type_double value) {
    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _columns; j++) {
            this->operator()(i, j) = value;
        }
    }
}

void Matrix::copy(Matrix& matrix) {
    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _columns; j++) {
            this->operator()(i, j) = matrix(i, j);
        }
    }
}

void Matrix::resize(size_t rows, size_t columns) {
    data.resize(rows*columns);
    _rows = rows;
    _columns = columns;
}