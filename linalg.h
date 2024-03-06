#ifndef LINALG_H
#define LINALG_H

#include <iostream>
#include <vector>
#include <initializer_list>
#include <array>
using namespace std;

class Vector {
    friend class Matrix;
private:
    std::vector<double> elements;
public:
    Vector(initializer_list<double> list);
    Vector(vector<double>& vec);
    size_t size() const;
    double getElementAt(size_t index) const;
    void setElementAt(size_t index, double value);
    Vector operator+(const Vector& other) const;
    Vector operator*(double scalar) const;
    double dot(const Vector& other) const;
    void print() const;
};

class Matrix {
private:
    vector<Vector> rows;

public: 
    //Constructors
    Matrix(initializer_list<Vector> list);
    Matrix(vector<Vector>& list);
    Matrix(size_t size);
    //Getters
    size_t numRows() const;
    size_t numCols() const;
    bool square() const;
    Vector getVector(size_t v) const;
    //setters
    void setElementAt(size_t row, size_t col, double value);
    //Basic Operators
    Matrix operator+(const Matrix& other) const;
    Matrix operator*(double scalar) const;
    Matrix operator*(const Matrix& other) const;
    //Elementary Row Operations
    void swapRows(size_t row1, size_t row2);
    void multiplyRow(size_t row, double scalar);
    void addRows(size_t targetRow, size_t sourceRow, double scalar);
    //More complicated matrix operations
    void rref();
    //Getters for data relevant to matrix
    void luDecomposition(Matrix& L, Matrix& U) const;
    double determinant() const;
    double trace() const;
    vector<Vector> rowSpace() const;
    vector<Vector> colSpace() const;
    vector<Vector> nullSpace() const;
    //Prints a representation of the matrix
    void print() const;
};

#endif // VECTORMATRIX_H
