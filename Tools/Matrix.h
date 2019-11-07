#ifndef CPP_RZ_PIC_MATRIX_H
#define CPP_RZ_PIC_MATRIX_H

#include <iostream>
#include <vector>
using namespace std;
//using matrix_type = double; // c++14 and higher
typedef double type_double; // c++0x


class Matrix {
private:
    size_t _rows;
    size_t _columns;
    //unique_ptr<matrix_type[]> data; // c++14 and higher
    vector<type_double> data;
public:
    Matrix() : _rows(0), _columns(0) {};
    Matrix(size_t rows, size_t columns);
    size_t rows() const;
    size_t columns() const;
    type_double* data_ptr();
    type_double& operator()(size_t row, size_t column);
    void print();
    void fill(const type_double value);
    void copy(Matrix& matrix);
    void resize(size_t rows, size_t columns);
};


#endif //CPP_RZ_PIC_MATRIX_H
