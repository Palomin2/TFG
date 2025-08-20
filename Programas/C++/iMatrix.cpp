/******************************************************************************* 
 *AUTHOR: Carlos Palomera Oliva
 *DATE: MARCH 2025
 *******************************************************************************/
#include <iostream>
#include <iomanip>
#include <vector>

#ifndef IMATRIX_H
#define IMATRIX_H

// This file is part of the qbLinAlg linear algebra library.

/*
MIT License
Copyright (c) 2023 Michael Bennett	

Permission is hereby granted, free of charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), to deal in the Software without restriction, 
including without limitation the rights to use, copy, modify, merge, publish, distribute, 
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is 
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or 
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING 
BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/




template <class T>
class iMatrix{
	public:
		// Define the various constructors.
    iMatrix();
    iMatrix(int nRows, int nCols);
    iMatrix(int nRows, int nCols, const T *inputData);
    iMatrix(const iMatrix<T> &inputMatrix);

    // And the destructor.
    ~iMatrix();

    // Configuration methods.
    bool Resize(int numRows, int numCols);
    void SetToIdentity();

    // Element access methods.
    T GetElement(int row, int col) const;
    std::vector<T> GetRow(int row) const;
    std::vector<T> GetCol(int col) const;    
    bool SetElement(int row, int col, T elementValue);
	bool SetCol(int col, std::vector<T> elements);
	bool SetRow(int row, std::vector<T> elements);
    int GetNumRows() const;
    int GetNumCols() const;
	iMatrix<T> GetSubMat(int rowIni, int rowEnd, int colIni, int colEnd) const;


    // Overload the assignment operator.
	iMatrix<T> operator= (const iMatrix<T> &rhs);


	// Write
    inline T& operator()(int i, int j) noexcept {
        return m_matrixData[i * m_nCols + j];
    };

    // Read
    inline const T& operator()(int i, int j) const noexcept {
        return m_matrixData[i * m_nCols + j];
    };


	iMatrix<T> invert() const { //Uses LU factorization with partial pivoting
        if (m_nRows != m_nCols) {
            throw std::runtime_error("Matrix must be square for inversion.");
        }

        int n = m_nRows;
        iMatrix<T> A(*this); 
        std::vector<int> piv(n);
        for (int i = 0; i < n; ++i) piv[i] = i;

        for (int k = 0; k < n; ++k) {
            T maxVal = std::abs(A(k, k));
            int pivotRow = k;
            for (int i = k + 1; i < n; ++i) {
                T val = std::abs(A(i, k));
                if (val > maxVal) {
                    maxVal = val;
                    pivotRow = i;
                }
            }
            if (maxVal == T(0)) {
                throw std::runtime_error("Matrix is singular.");
            }

            if (pivotRow != k) {
                for (int j = 0; j < n; ++j) {
                    std::swap(A(k, j), A(pivotRow, j));
                }
                std::swap(piv[k], piv[pivotRow]);
            }

            for (int i = k + 1; i < n; ++i) {
                A(i, k) /= A(k, k);
                T factor = A(i, k);
                for (int j = k + 1; j < n; ++j) {
                    A(i, j) -= factor * A(k, j);
                }
            }
        }

        iMatrix<T> inv(n, n, T(0));
        for (int col = 0; col < n; ++col) {
            std::vector<T> b(n, T(0));
            b[piv[col]] = T(1);

            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < i; ++j) {
                    b[i] -= A(i, j) * b[j];
                }
            }

            for (int i = n - 1; i >= 0; --i) {
                for (int j = i + 1; j < n; ++j) {
                    b[i] -= A(i, j) * b[j];
                }
                b[i] /= A(i, i);
            }

            for (int i = 0; i < n; ++i) {
                inv(i, col) = b[i];
            }
        }

        return inv;
    }
    // Overload +, - and * operators (friends).
    template <class U> friend iMatrix<U> operator+ (const iMatrix<U>& lhs, const iMatrix<U>& rhs);
    template <class U> friend iMatrix<U> operator+ (const U& lhs, const iMatrix<U>& rhs);
    template <class U> friend iMatrix<U> operator+ (const iMatrix<U>& lhs, const U& rhs);
        
    template <class U> friend iMatrix<U> operator- (const iMatrix<U>& lhs, const iMatrix<U>& rhs);
    template <class U> friend iMatrix<U> operator- (const U& lhs, const iMatrix<U>& rhs);
    template <class U> friend iMatrix<U> operator- (const iMatrix<U>& lhs, const U& rhs);
        
    template <class U> friend iMatrix<U> operator* (const iMatrix<U>& lhs, const iMatrix<U>& rhs);
    template <class U> friend iMatrix<U> operator* (const U& lhs, const iMatrix<U>& rhs);
    template <class U> friend iMatrix<U> operator* (const iMatrix<U>& lhs, const U& rhs);

	

    void PrintMatrix();
    void PrintMatrix(int precision);

    bool IsSquare();
	private:
		int Sub2Ind(int row, int col) const;
    private:
    T *m_matrixData;
    int m_nRows, m_nCols, m_nElements;
};   


/* **************************************************************************************************
CONSTRUCTOR / DESTRUCTOR FUNCTIONS
/* *************************************************************************************************/
// The default constructor.
template <class T>
iMatrix<T>::iMatrix()
{
  m_nRows = 1;
  m_nCols = 1;
  m_nElements = 1;
  m_matrixData = nullptr;
}

// Construct empty matrix (all elements 0)
template <class T>
iMatrix<T>::iMatrix(int nRows, int nCols)
{
  m_nRows = nRows;
  m_nCols = nCols;
  m_nElements = m_nRows * m_nCols;
  m_matrixData = new T[m_nElements];
  for (int i = 0; i < m_nElements; i++) {
	if constexpr (std::is_arithmetic_v<T>) {
		// Si T es un tipo numérico (int, float, double, etc.)
		m_matrixData[i] = 0.0f;
	} else {
		// Si T es un contenedor (std::vector<float>, etc.)
		m_matrixData[i] = T(); // Constructor por defecto (vector vacío)
	}
}
}

// Construct from const linear array.
template <class T>
iMatrix<T>::iMatrix(int nRows, int nCols, const T *inputData)
{
	m_nRows = nRows;
	m_nCols = nCols;
	m_nElements = m_nRows * m_nCols;
	m_matrixData = new T[m_nElements];
	for (int i=0; i<m_nElements; i++)
	  m_matrixData[i] = inputData[i];
}

// The copy constructor.
template <class T>
iMatrix<T>::iMatrix(const iMatrix<T> &inputMatrix)
{
	m_nRows = inputMatrix.m_nRows;
	m_nCols = inputMatrix.m_nCols;
	m_nElements = inputMatrix.m_nElements;
	
	m_matrixData = new T[m_nElements];
	for (int i=0; i<m_nElements; i++)
		m_matrixData[i] = inputMatrix.m_matrixData[i];
}

template <class T>
iMatrix<T>::~iMatrix()
{
	// Destructor.
	if (m_matrixData)
		delete[] m_matrixData;
	
	m_matrixData = nullptr;
}
/* **************************************************************************************************
CONFIGURATION FUNCTIONS
/* *************************************************************************************************/
template <class T>
bool iMatrix<T>::Resize(int numRows, int numCols)
{
	m_nRows = numRows;
	m_nCols = numCols;
	m_nElements = (m_nRows * m_nCols);
	delete[] m_matrixData;
	m_matrixData = new T[m_nElements];
	if (m_matrixData != nullptr)
	{
		for (int i=0; i<m_nElements; i++)
			m_matrixData[i] = 0.0f;

		return true;
	}
	else
	{
		return false;
	}
}
template <class T>
void iMatrix<T>::SetToIdentity()
{
	if (!IsSquare())
		throw std::invalid_argument("Cannot form an identity matrix that is not square.");
		
	for (int row=0; row<m_nRows; ++row)
	{
		for (int col=0; col<m_nCols; ++col)
		{
			if (col == row)
				m_matrixData[Sub2Ind(row,col)] = 1.0;
			else
				m_matrixData[Sub2Ind(row,col)] = 0.0;
		}
	}
}

template <class T>
int iMatrix<T>::GetNumRows() const
{
	return m_nRows;
}

template <class T>
int iMatrix<T>::GetNumCols() const
{
	return m_nCols;
}
/* **************************************************************************************************
ELEMENT FUNCTIONS
/* *************************************************************************************************/
template <class T>
T iMatrix<T>::GetElement(int row, int col) const
{
	int linearIndex = Sub2Ind(row, col);
	if (linearIndex >= 0)
		return m_matrixData[linearIndex];
	else{
        if constexpr (std::is_arithmetic_v<T>) {
            return 0; // Para tipos numéricos (int, float, double)
        } else {
            return T(); // Para contenedores (std::vector<float>, etc.)
        }
    }
}

template <class T>
std::vector<T> iMatrix<T>::GetRow(int row) const
{   
    int linearIndex = Sub2Ind(row, 0);
    std::vector<T> AuxRow;
    for(unsigned int i=0;i<this->GetNumCols();i++){
        AuxRow.push_back(m_matrixData[linearIndex]);
        linearIndex++;
    }
    return AuxRow;

}

template <class T>
bool iMatrix<T>::SetElement(int row, int col, T elementValue)
{
	int linearIndex = Sub2Ind(row, col);
	if (linearIndex >= 0)
	{
		m_matrixData[linearIndex] = elementValue;
		return true;
	} 
	else 
	{
		return false;
	}
}

template <class T>
bool iMatrix<T>::SetCol(int col, std::vector<T> elements){
	int linearIndex = Sub2Ind(0, col);
	int aux =this->GetNumRows();
	if (linearIndex >= 0){
		for(unsigned int i=0;i<aux;i++){
			(*this)(i,col)=elements[i];
			//this->SetElement(i,col,elements[i]);
			linearIndex=linearIndex+aux;
		}
		return true;
	}
	else{
		return false;
	}
}

template <class T>
bool iMatrix<T>::SetRow(int row, std::vector<T> elements){
	int linearIndex = Sub2Ind(row, 0);
	int aux =this->GetNumCols();
	if (linearIndex >= 0){
		for(unsigned int i=0;i<aux;i++){
			(*this)(row, i) = elements[i];
			//this->SetElement(row,i,elements[i]);
			linearIndex=linearIndex+aux;
		}
		return true;
	}
	else{
		return false;
	}

}

template <class T>
std::vector<T> iMatrix<T>::GetCol(int col) const
{   
    int linearIndex = Sub2Ind(0, col);
    std::vector<T> AuxCol;
    int aux =this->GetNumCols();
	int maxiter = this->GetNumRows();
    for(unsigned int i=0; i<maxiter; i++){
        AuxCol.push_back(m_matrixData[linearIndex]);
        linearIndex=linearIndex+aux;
    }
    return AuxCol;

}

template <class T>
iMatrix<T> iMatrix<T>::GetSubMat(int rowIni, int rowEnd, int colIni, int colEnd) const
{   
	int numRows= rowEnd-rowIni+1;
	int numCols= colEnd-colIni+1;

	int aux =this->GetNumCols();
    iMatrix<T> AuxMat(numRows, numCols);

	for(unsigned int i=0; i< numRows; i++){
		for(unsigned int j=0; j<numCols; j++){
			AuxMat(i,j)=(*this)(i+rowIni, j+colIni);
		}
	}

	return AuxMat;
}

/* **************************************************************************************************
THE + OPERATOR
/* *************************************************************************************************/
// matrix + matrx.
template <class T>
iMatrix<T> operator+ (const iMatrix<T>& lhs, const iMatrix<T>& rhs)
{
	int numRows = lhs.m_nRows;
	int numCols = lhs.m_nCols;
	int numElements = numRows * numCols;
	T *tempResult = new T[numElements];
	for (int i=0; i<numElements; i++)
		tempResult[i] = lhs.m_matrixData[i] + rhs.m_matrixData[i];
		
	iMatrix<T> result(numRows, numCols, tempResult);
	delete[] tempResult;
	return result;
}

// scaler + matrix
template <class T>
iMatrix<T> operator+ (const T& lhs, const iMatrix<T>& rhs)
{
	int numRows = rhs.m_nRows;
	int numCols = rhs.m_nCols;
	int numElements = numRows * numCols;
	T *tempResult = new T[numElements];
	for (int i=0; i<numElements; ++i)
		tempResult[i] = lhs + rhs.m_matrixData[i];
		
	iMatrix<T> result(numRows, numCols, tempResult);
	delete[] tempResult;
	return result;
}

// matrix + scaler
template <class T>
iMatrix<T> operator+ (const iMatrix<T>& lhs, const T& rhs)
{
	int numRows = lhs.m_nRows;
	int numCols = lhs.m_nCols;
	int numElements = numRows * numCols;
	T *tempResult = new T[numElements];
	for (int i=0; i<numElements; ++i)
		tempResult[i] = lhs.m_matrixData[i] + rhs;
		
	iMatrix<T> result(numRows, numCols, tempResult);
	delete[] tempResult;
	return result;
}

/* **************************************************************************************************
THE - OPERATOR
/* *************************************************************************************************/
// matrix - matrix
template <class T>
iMatrix<T> operator- (const iMatrix<T>& lhs, const iMatrix<T>& rhs)
{
	int numRows = lhs.m_nRows;
	int numCols = lhs.m_nCols;
	int numElements = numRows * numCols;
	T *tempResult = new T[numElements];
	for (int i=0; i<numElements; i++)
		tempResult[i] = lhs.m_matrixData[i] - rhs.m_matrixData[i];
		
	iMatrix result(numRows, numCols, tempResult);
	delete[] tempResult;
	return result;    
}

// scaler - matrix
template <class T>
iMatrix<T> operator- (const T& lhs, const iMatrix<T>& rhs)
{
	int numRows = rhs.m_nRows;
	int numCols = rhs.m_nCols;
	int numElements = numRows * numCols;
	T *tempResult = new T[numElements];
	for (int i=0; i<numElements; ++i)
		tempResult[i] = lhs - rhs.m_matrixData[i];
		
	iMatrix<T> result(numRows, numCols, tempResult);
	delete[] tempResult;
	return result;
}

// matrix - scaler
template <class T>
iMatrix<T> operator- (const iMatrix<T>& lhs, const T& rhs)
{
	int numRows = lhs.m_nRows;
	int numCols = lhs.m_nCols;
	int numElements = numRows * numCols;
	T *tempResult = new T[numElements];
	for (int i=0; i<numElements; ++i)
		tempResult[i] = lhs.m_matrixData[i] - rhs;
		
	iMatrix<T> result(numRows, numCols, tempResult);
	delete[] tempResult;
	return result;
}

// scaler * matrix
template <class T>
iMatrix<T> operator* (const T& lhs, const iMatrix<T>& rhs)
{
	int numRows = rhs.m_nRows;
	int numCols = rhs.m_nCols;
	int numElements = numRows * numCols;
	T *tempResult = new T[numElements];
	for (int i=0; i<numElements; ++i)
		tempResult[i] = lhs * rhs.m_matrixData[i];
		
	iMatrix<T> result(numRows, numCols, tempResult);
	delete[] tempResult;
	return result;
}

// matrix * scaler
template <class T>
iMatrix<T> operator* (const iMatrix<T>& lhs, const T& rhs)
{
	int numRows = lhs.m_nRows;
	int numCols = lhs.m_nCols;
	int numElements = numRows * numCols;
	T *tempResult = new T[numElements];
	for (int i=0; i<numElements; ++i)
		tempResult[i] = lhs.m_matrixData[i] * rhs;
		
	iMatrix<T> result(numRows, numCols, tempResult);
	delete[] tempResult;
	return result;
}

// matrix * matrix
template <class T>
iMatrix<T> operator* (const iMatrix<T>& lhs, const iMatrix<T>& rhs)
{
	int r_numRows = rhs.m_nRows;
	int r_numCols = rhs.m_nCols;
	int l_numRows = lhs.m_nRows;
	int l_numCols = lhs.m_nCols;

	if (l_numCols == r_numRows)
	{
		// This is the standard matrix multiplication condition.
		// The output will be the same size as the RHS.
		T *tempResult = new T[lhs.m_nRows * rhs.m_nCols];

		// Loop through each row of the LHS.
		for (int lhsRow=0; lhsRow<l_numRows; lhsRow++)
		{
			// Loop through each column on the RHS.
			for (int rhsCol=0; rhsCol<r_numCols; rhsCol++)
			{
				T elementResult = static_cast<T>(0.0);
				// Loop through each element of this LHS row.
				for (int lhsCol=0; lhsCol<l_numCols; lhsCol++)
				{
					// Compute the LHS linear index.
					int lhsLinearIndex = (lhsRow * l_numCols) + lhsCol;

					// Compute the RHS linear index (based on LHS col).
					// rhs row number equal to lhs column number.
					int rhsLinearIndex = (lhsCol * r_numCols) + rhsCol;

					// Perform the calculation on these elements.
					elementResult += (lhs.m_matrixData[lhsLinearIndex] * rhs.m_matrixData[rhsLinearIndex]);
				}

				// Store the result.
				int resultLinearIndex = (lhsRow * r_numCols) + rhsCol;
				tempResult[resultLinearIndex] = elementResult;
			}		
		}
		iMatrix<T> result(l_numRows, r_numCols, tempResult);
		delete[] tempResult;
		return result;
	}
	else
	{
		iMatrix<T> result(1, 1);
		return result;
	}
}

template <class T>
iMatrix<T> iMatrix<T>::operator= (const iMatrix<T> &rhs)
{
	// Make sure we're not assigning to ourself.
	if (this != &rhs)
	{
		// If the dimensions are the same, we only need to copy the elements,
		//	there is no need to delete and re-allocate memory.
		if ((m_nRows == rhs.m_nRows) && (m_nCols == rhs.m_nCols))
		{
			for (int i=0; i<m_nElements; ++i)
				m_matrixData[i] = rhs.m_matrixData[i];
		}
		else
		{
			m_nRows = rhs.m_nRows;
			m_nCols = rhs.m_nCols;
			m_nElements = rhs.m_nElements;
			
			if (m_matrixData)
				delete[] m_matrixData;
			
			m_matrixData = new T[m_nElements];
			for (int i=0; i<m_nElements; i++)
				m_matrixData[i] = rhs.m_matrixData[i];	
		}
	}
	
	return *this;
}

// A simple function to print a matrix to stdout.
template <class T>
void iMatrix<T>::PrintMatrix()
{
	int nRows = this->GetNumRows();
	int nCols = this->GetNumCols();
	for (int row = 0; row<nRows; ++row)
  {
	  for (int col = 0; col<nCols; ++col)
    {
	    std::cout << std::fixed << std::setprecision(3) << (*this)(row, col) << "  ";
    }
	std::cout << std::endl;
	}    
}

// A simple function to print a matrix to stdout, with specified precision.
template <class T>
void iMatrix<T>::PrintMatrix(int precision)
{
	int nRows = this->GetNumRows();
	int nCols = this->GetNumCols();
	for (int row = 0; row<nRows; ++row)
  {
	  for (int col = 0; col<nCols; ++col)
    {
	    std::cout << std::fixed << std::setprecision(precision) << (*this)(row, col) << "  ";
    }
	std::cout << std::endl;
	}    
}



/* **************************************************************************************************
PRIVATE FUNCTIONS
/* *************************************************************************************************/
// Function to return the linear index corresponding to the supplied row and column values.
template <class T>
int iMatrix<T>::Sub2Ind(int row, int col) const
{
	if ((row < m_nRows) && (row >= 0) && (col < m_nCols) && (col >= 0))
		return (row * m_nCols) + col;
	else
		return -1;
}

// Function to test whether the matrix is square.
template <class T>
bool iMatrix<T>::IsSquare()
{
	if (m_nCols == m_nRows)
		return true;
	else
		return false;
}

#endif