#include "linalg.h"

/// @brief Constructs a Vector
/// @param list A set of scalars
Vector::Vector(std::initializer_list<double> list) : elements(list) {}

/// @brief Constructs a vector from a set of scalars (for function returns)
/// @param vec a set of scalars
Vector::Vector(vector<double> &vec) : elements(vec) {}

/// @brief Getter for size of vector
/// @return size_t size
size_t Vector::size() const
{
    return elements.size();
}

/// @brief Gets the element at the specified index
/// @param index target index
/// @return The element at the specified index
double Vector::getElementAt(size_t index) const
{
    if (index > elements.size())
    {
        throw invalid_argument("Index bigger than vector!");
    }
    return this->elements[index];
}

/// @brief Sets the element at the specified index in a vector
/// @param index The index to modify
/// @param value The new value to set
void Vector::setElementAt(size_t index, double value)
{
    if (index >= elements.size())
    {
        throw invalid_argument("Index bigger than vector!");
    }
    elements[index] = value;
}

/// @brief Vector addition operator
/// @param other The vector being added
/// @return a new vector representing the sum of the two vectors
Vector Vector::operator+(const Vector &other) const
{
    if (elements.size() != other.size())
    {
        throw invalid_argument("Vector size mismatch!");
    }

    Vector result({});
    for (size_t i = 0; i < elements.size(); ++i)
    {
        result.elements.push_back(this->elements[i] + other.elements[i]);
    }
    return result;
}
/// @brief Performs scalar multiplication on the given vector
/// @param scalar The scalar for product
/// @return New vector representing the product of the vector and the scalar
Vector Vector::operator*(double scalar) const
{
    Vector result({});
    for (double elem : elements)
    {
        result.elements.push_back(elem * scalar);
    }
    return result;
}

/// @brief Calculates the dot product of two vectors
/// @param other The second vector
/// @return the resulting dot product
double Vector::dot(const Vector &other) const
{
    if (this->size() != other.size())
    {
        throw invalid_argument("Vector size mismatch!");
    }

    double result = 0;
    for (size_t i = 0; i < this->size(); i++)
    {
        result += this->elements[i] * other.elements[i];
    }
    return result;
}

/// @brief Prints a representation of the vector
void Vector::print() const
{
    cout << "(";
    for (size_t i = 0; i < this->size(); i++)
    {
        cout << this->elements[i];
        if (i < this->size() - 1)
        {
            cout << ", ";
        }
    }
    cout <<")"<<endl;
}

/// @brief Constructs a matrix from a set of vectors
/// @param list a set of vectors
Matrix::Matrix(initializer_list<Vector> list) : rows(list) {}

/// @brief Constructs a matrix from a reference to a set of vectors (for building matrix for function returns)
/// @param list The set of vectors
Matrix::Matrix(vector<Vector> &list) : rows(list) {}

/// @brief Creates identity matrix of specified size
/// @param size The size of the identity matrix
Matrix::Matrix(size_t size)
{
    if (size == 0)
    {
        throw invalid_argument("Cannot create matrix of size 0");
    }
    std::vector<double> row(size, 0);
    for (size_t i = 0; i < size; ++i)
    {
        row[i] = 1;
        rows.push_back(Vector(row));
        row[i] = 0;
    }
}

/// @brief Gets the number of rows of the matrix
/// @return number of rows
size_t Matrix::numRows() const
{
    return rows.size();
}

/// @brief Gets the number of columns of a matrix
/// @return
size_t Matrix::numCols() const
{
    return rows.empty() ? 0 : rows[0].size();
}

/// @brief Returns whether or not the matrix is square
/// @return boolean flag
bool Matrix::square() const
{
    return this->numRows() == this->numCols();
}

/// @brief Gets a vector representing a row of the matrix
/// @param index The index of the row
/// @return The vector at the specified row
Vector Matrix::getVector(size_t index) const
{
    if (index >= numRows())
    {
        throw invalid_argument("Index out of bounds.");
    }
    return rows[index];
}

/// @brief Sets the element at the specified row and column in a matrix
/// @param row The row of the element
/// @param col The column of the element
/// @param value The new value to set
void Matrix::setElementAt(size_t row, size_t col, double value)
{
    if (row >= numRows() || col >= numCols())
    {
        throw invalid_argument("Index out of range.");
    }
    rows[row].setElementAt(col, value);
}

/// @brief Performs addition between two matricies
/// @param other The second matrix
/// @return The sum of the two matricies
Matrix Matrix::operator+(const Matrix &other) const
{
    if (this->numRows() != other.numRows() || this->numCols() != other.numCols())
    {
        throw invalid_argument("Matrix Size Mismatch!");
    }

    vector<Vector> resultRows;
    for (size_t i = 0; i < this->numRows(); ++i)
    {
        resultRows.push_back(this->rows[i] + other.rows[i]);
    }
    return Matrix(resultRows);
}

/// @brief Multiplies the matrix by a given scalar
/// @param scalar a scalar
/// @return The product matrix
Matrix Matrix::operator*(double scalar) const
{
    std::vector<Vector> resultRows;
    for (const Vector &row : rows)
    {
        resultRows.push_back(row * scalar);
    }
    return Matrix(resultRows);
}

/// @brief Multiplies two matricies
/// @param other The second matrix
/// @return A new matrix representing the product of the two matricies
Matrix Matrix::operator*(const Matrix &other) const
{
    if (this->numCols() != other.numRows())
    {
        throw invalid_argument("Matrix Size Mismatch!");
    }

    vector<Vector> resultRows;
    for (size_t i = 0; i < this->numRows(); ++i)
    {
        vector<double> newRowElements;
        for (size_t j = 0; j < other.numCols(); ++j)
        {
            double sum = 0;
            for (size_t k = 0; k < this->numCols(); ++k)
            {
                sum += this->rows[i].getElementAt(k) * other.rows[k].getElementAt(j);
            }
            newRowElements.push_back(sum);
        }
        resultRows.push_back(Vector(newRowElements));
    }
    return Matrix(resultRows);
}

/// @brief The first elementary row operation: Swapping two rows
/// @param row1 The index of the first row
/// @param row2 The index of the second row
void Matrix::swapRows(size_t row1, size_t row2)
{
    if (row1 >= numRows() || row2 >= numRows())
    {
        throw invalid_argument("Row index out of bounds");
    }
    std::swap(rows[row1], rows[row2]);
}

/// @brief The second elementary row operation: Multiplies a row by a scalar
/// @param row The target row
/// @param scalar some scalar
void Matrix::multiplyRow(size_t row, double scalar)
{
    if (row >= numRows())
    {
        throw invalid_argument("Row index out of bounds");
    }
    if (scalar == 0)
    {
        throw invalid_argument("Scalar cannot be zero");
    }
    rows[row] = rows[row] * scalar;
}

/// @brief The third elementary row operation: Adding a multiple of a row to another
/// @param targetRow The target row
/// @param sourceRow The original row
/// @param scalar The scalar to multiply the original by
void Matrix::addRows(size_t targetRow, size_t sourceRow, double scalar)
{
    if (targetRow >= numRows() || sourceRow >= numRows())
    {
        throw invalid_argument("Row index out of bounds");
    }

    if (rows[targetRow].size() != rows[sourceRow].size())
    {
        throw invalid_argument("Row size mismatch");
    }

    for (size_t i = 0; i < rows[targetRow].size(); ++i)
    {
        double newValue = rows[targetRow].getElementAt(i) + rows[sourceRow].getElementAt(i) * scalar;
        rows[targetRow].elements[i] = newValue;
    }
}

/// @brief Places the given matrix in Reduced Row Echelon Form
void Matrix::rref()
{
    size_t numRows = this->numRows();
    size_t numCols = this->numCols();
    size_t lead = 0;

    for (size_t r = 0; r < numRows; ++r)
    {
        if (lead >= numCols)
        {
            break;
        }
        size_t i = r;
        while (rows[i].getElementAt(lead) == 0)
        {
            ++i;
            if (i == numRows)
            {
                i = r;
                ++lead;
                if (lead == numCols)
                {
                    return;
                }
            }
        }
        swapRows(i, r);

        double div = rows[r].getElementAt(lead);
        if (div != 0)
        {
            multiplyRow(r, 1 / div);
        }

        for (size_t j = 0; j < numRows; ++j)
        {
            if (j != r)
            {
                double mult = -rows[j].getElementAt(lead);
                addRows(j, r, mult);
            }
        }
        ++lead;
    }
}

/// @brief Calculates the trace of the matrix
/// @return The trace of the matrix
double Matrix::trace() const
{
    if (!this->square())
    {
        throw invalid_argument("Trace is only defined for square matrices.");
    }

    double trace = 0;
    for (size_t i = 0; i < this->numRows(); ++i)
    {
        trace += this->rows[i].getElementAt(i);
    }

    return trace;
}

/// @brief Returns a basis for the row space of this matrix
/// @return a span of vectors
vector<Vector> Matrix::rowSpace() const
{
    Matrix rrefMatrix = *this;
    rrefMatrix.rref();

    vector<Vector> basis;
    for (const auto &row : rrefMatrix.rows)
    {
        bool nonZeroRow = false;
        for (double elem : row.elements)
        {
            if (elem != 0)
            {
                nonZeroRow = true;
                break;
            }
        }
        if (nonZeroRow)
        {
            basis.push_back(row);
        }
    }
    return basis;
}

/// @brief Returns a basis for the column space of this matrix
/// @return a span of vectors
vector<Vector> Matrix::colSpace() const
{
    Matrix rrefMatrix = *this;
    rrefMatrix.rref(); // Reduce to RREF

    vector<Vector> basis;
    for (size_t i = 0; i < rrefMatrix.numCols(); ++i)
    {
        bool isLeadingOne = false;
        for (size_t j = 0; j < rrefMatrix.numRows() && !isLeadingOne; ++j)
        {
            if (rrefMatrix.rows[j].getElementAt(i) == 1)
            {
                bool isOnlyNonZero = true;
                for (size_t k = 0; k < rrefMatrix.numRows(); ++k)
                {
                    if (k != j && rrefMatrix.rows[k].getElementAt(i) != 0)
                    {
                        isOnlyNonZero = false;
                        break;
                    }
                }
                if (isOnlyNonZero)
                    isLeadingOne = true;
            }
        }
        if (isLeadingOne)
        {
            vector<double> colVec(this->numRows());
            for (size_t j = 0; j < this->numRows(); ++j)
            {
                colVec[j] = this->rows[j].getElementAt(i);
            }
            basis.push_back(Vector(colVec));
        }
    }
    return basis;
}

/// @brief Returns a basis for the nullspace space of this matrix
/// @return a span of vectors
vector<Vector> Matrix::nullSpace() const
{
    Matrix rrefMatrix = *this;
    rrefMatrix.rref();

    size_t n = rrefMatrix.numCols();
    vector<Vector> basis;

    vector<bool> isFreeVariable(n, true);
    for (const auto &row : rrefMatrix.rows)
    {
        for (size_t i = 0; i < n; ++i)
        {
            if (row.getElementAt(i) == 1)
            {
                isFreeVariable[i] = false;
                break;
            }
        }
    }

    for (size_t i = 0; i < n; ++i)
    {
        if (isFreeVariable[i])
        {
            std::vector<double> nullSpaceVector(n, 0);
            nullSpaceVector[i] = 1;

            for (const auto &row : rrefMatrix.rows)
            {
                for (size_t j = 0; j < n; ++j)
                {
                    if (row.getElementAt(j) == 1)
                    {

                        double sum = 0;
                        for (size_t k = j + 1; k < n; ++k)
                        {
                            sum += row.getElementAt(k) * nullSpaceVector[k];
                        }
                        nullSpaceVector[j] = -sum;
                        break;
                    }
                }
            }
            basis.push_back(Vector(nullSpaceVector));
        }
    }
    return basis;
}

/// @brief Performs imperfect LU decomposition TODO: Partial pivot implementation
/// @param L References to lower matrix to be created
/// @param U Reference to upper matrix to be created
void Matrix::luDecomposition(Matrix &L, Matrix &U) const
{
    if (!this->square())
    {
        throw invalid_argument("Matrix must be square.");
    }
    size_t n = this->numRows();
    L = Matrix(n);
    U = *this;

    for (size_t i = 0; i < n; i++)
    {
        for (size_t k = i + 1; k < n; k++)
        {
            if (U.rows[i].getElementAt(i) == 0)
            {
                throw invalid_argument("Zero pivot encountered.");
            }

            double m = U.rows[k].getElementAt(i) / U.rows[i].getElementAt(i);
            L.rows[k].setElementAt(i, m);
            for (size_t j = 0; j < n; j++)
            {
                U.rows[k].setElementAt(j, U.rows[k].getElementAt(j) - m * U.rows[i].getElementAt(j));
            }
        }
    }
    for (size_t i = 0; i < n; i++)
    {
        L.rows[i].setElementAt(i, 1);
    }
}

/// @brief Calculates the determinant of this matrix, given that it is square
/// @return Determinant of this matrix
double Matrix::determinant() const
{
    if (!this->square())
    {
        throw invalid_argument("Determinant is only defined for square matrices.");
    }

    size_t n = this->numRows();
    Matrix L = Matrix(n);
    Matrix U = Matrix(n);
    this->luDecomposition(L, U);
    double det = 1;
    for (size_t i = 0; i < this->numRows(); i++)
    {
        det *= U.rows[i].getElementAt(i);
    }
    return det;
}

/// @brief Prints a representation of this matrix
void Matrix::print() const {
    for (const auto& row : rows) {
        std::cout << "(";
        for (size_t i = 0; i < row.size(); ++i) {
            std::cout << row.getElementAt(i);
            if (i < row.size() - 1) std::cout << ", ";
        }
        std::cout << ")" << std::endl;
    }
}