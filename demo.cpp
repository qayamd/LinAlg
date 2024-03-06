#include <iostream>
#include <vector>
#include "linalg.h"

using namespace std;

int main() {
    // Demonstrate vector operations
    cout << "Vector Operations:" << endl;
    Vector v1({1, 2, 3});
    Vector v2({4, 5, 6});
    Vector v3 = v1 + v2; // Vector addition
    cout << "v1 + v2 = ";
    v3.print();
    cout << endl;

    double scalar = 2.0;
    Vector v4 = v1 * scalar; // Scalar multiplication
    cout << "v1 * " << scalar << " = ";
    v4.print();
    cout << endl;

    double dotProduct = v1.dot(v2); // Dot product
    cout << "Dot product of v1 and v2 = " << dotProduct << endl << endl;

    // Demonstrate matrix operations
    cout << "Matrix Operations:" << endl;
    Matrix m1({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
    Matrix m2({{9, 8, 7}, {6, 5, 4}, {3, 2, 1}});
    Matrix m3 = m1 + m2; // Matrix addition
    cout << "m1 + m2 = " << endl;
    m3.print();
    cout << endl;

    Matrix m4 = m1 * scalar; // Scalar multiplication
    cout << "m1 * " << scalar << " = " << endl;
    m4.print();
    cout << endl;

    Matrix m5 = m1 * m2; // Matrix multiplication
    cout << "m1 * m2 = " << endl;
    m5.print();
    cout << endl;

    double det = m1.determinant(); // Determinant
    cout << "Determinant of m1 = " << det << endl << endl;

    // Row space, column space, and null space
    cout << "Row space of m1:" << endl;
    vector<Vector> rowSpace = m1.rowSpace();
    for (auto& vec : rowSpace) {
        vec.print();
        cout << endl;
    }

    cout << "Column space of m1:" << endl;
    vector<Vector> colSpace = m1.colSpace();
    for (auto& vec : colSpace) {
        vec.print();
        cout << endl;
    }

    cout << "Null space of m1:" << endl;
    vector<Vector> nullSpace = m1.nullSpace();
    for (auto& vec : nullSpace) {
        vec.print();
        cout << endl;
    }

    return 0;
}
