/******************************************************************************* 
 *AUTHOR: Carlos Palomera Oliva
 *DATE: MARCH 2025
 *******************************************************************************/
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include "C:\Librerias\eigen-3.4.0\Eigen\Dense"   //windows version
//#include "/home/carlo/Librerias/eigen-3.4.0/Eigen/Dense"  //linux version
//#include "/home/carlo/Librerias/eigen-3.4.0/Eigen/Sparse" //linux version
#include <iomanip>
#include <vector>
#include <cmath>
#include <numeric>
#include <type_traits>

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
	iMatrix(int nRows, int nCols, const Eigen::VectorXd &inputData);
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

        iMatrix<T> inv(n, n);
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

	template <class U> friend iMatrix<U> operator* (const std::vector<U>& v1, const std::vector<U>& v2);

	iMatrix<T>& operator+= (const iMatrix<T>& rhs);
	iMatrix<T>& operator*= (const T& rhs);

    void PrintMatrix();
    void PrintMatrix(int precision);
	std::vector<T> Vectorize() const;
	iMatrix<T> Transpose() const;


    bool IsSquare();
	T Determinant() const;

	// Función estática para calcular el Producto Punto entre dos std::vector<T>.
	static T dot_product (const std::vector<T>& v1, const std::vector<T>& v2);

	private:
		int Sub2Ind(int row, int col) const;
    private:
    std::vector<T> m_matrixData;
    int m_nRows, m_nCols;
};   


/* **************************************************************************************************
CONSTRUCTOR / DESTRUCTOR FUNCTIONS
/* *************************************************************************************************/
template <class T>
iMatrix<T>::iMatrix() : m_nRows(1), m_nCols(1)
{
  // Inicializamos el vector con 1 elemento (por defecto 0 o T())
  m_matrixData.resize(1);
}

// Construct empty matrix (all elements 0)
template <class T>
iMatrix<T>::iMatrix(int nRows, int nCols) : m_nRows(nRows), m_nCols(nCols)
{
  int nElements = m_nRows * m_nCols;
  if constexpr (std::is_arithmetic_v<T>) {
		// Redimensionar e inicializar a 0.0f
		m_matrixData.resize(nElements, static_cast<T>(0.0f));
	} else {
		// Redimensionar e inicializar con el constructor por defecto T()
		m_matrixData.resize(nElements); 
	}
}
// Construct from Eigen::VectorXd
template <class T>
iMatrix<T>::iMatrix(int nRows, int nCols, const Eigen::VectorXd& inputData)
    : m_nRows(nRows), m_nCols(nCols)
{
    int nElements = m_nRows * m_nCols;
    if (inputData.size() != nElements)
        throw std::invalid_argument("Eigen::VectorXd no tiene el tamaño correcto para llenar la matriz.");

    m_matrixData.resize(nElements);
    for (int i = 0; i < nElements; ++i)
        m_matrixData[i] = static_cast<T>(inputData[i]);
}
// Construct from const linear array.
template <class T>
iMatrix<T>::iMatrix(int nRows, int nCols, const T *inputData) : m_nRows(nRows), m_nCols(nCols)
{
	int nElements = m_nRows * m_nCols;
	// Usamos el constructor de rango de std::vector para la copia eficiente.
	m_matrixData = std::vector<T>(inputData, inputData + nElements);
}

// The copy constructor.
template <class T>
iMatrix<T>::iMatrix(const iMatrix<T> &inputMatrix) 
    : m_nRows(inputMatrix.m_nRows), 
      m_nCols(inputMatrix.m_nCols), 
      m_matrixData(inputMatrix.m_matrixData) 
{
	// Cuerpo del constructor vacío (RAII).
}

template <class T>
iMatrix<T>::~iMatrix()
{
	// El destructor de std::vector se llama automáticamente (RAII).
}
/* **************************************************************************************************
CONFIGURATION FUNCTIONS
/* *************************************************************************************************/
template <class T>
bool iMatrix<T>::Resize(int numRows, int numCols)
{
	m_nRows = numRows;
	m_nCols = numCols;
	int nElements = (m_nRows * m_nCols);
	
	
	if constexpr (std::is_arithmetic_v<T>) {
		m_matrixData.resize(nElements, static_cast<T>(0.0f));
	} else {
		m_matrixData.resize(nElements); 
	}

	return true; 
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
    if (col < 0 || col >= m_nCols || static_cast<int>(elements.size()) != m_nRows)
        return false;
    for (int i = 0; i < m_nRows; ++i) {
        m_matrixData[i * m_nCols + col] = elements[i];
    }
    return true;
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
    std::vector<T> result(m_nRows);
    for (int i = 0; i < m_nRows; ++i) {
        result[i] = m_matrixData[i * m_nCols + col];
    }
    return result;
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
	int numElements = lhs.m_matrixData.size();
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
	int numElements = rhs.m_matrixData.size();
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
	int numElements = lhs.m_matrixData.size();
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
	int numElements = lhs.m_matrixData.size();
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
	int numElements = rhs.m_matrixData.size();
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
	int numElements = lhs.m_matrixData.size();
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
    iMatrix<T> result(rhs.m_nRows, rhs.m_nCols);
    std::transform(
        rhs.m_matrixData.begin(), 
        rhs.m_matrixData.end(),  
        result.m_matrixData.begin(),
        [&lhs](const T& element) { return lhs * element; } 
    );
    return result;
}

/* Old product operator for the scaler

template <class T>
iMatrix<T> operator* (const T& lhs, const iMatrix<T>& rhs)
{
	int numRows = rhs.m_nRows;
	int numCols = rhs.m_nCols;
	int numElements = rhs.m_matrixData.size();
	T *tempResult = new T[numElements];
	for (int i=0; i<numElements; ++i)
		tempResult[i] = lhs * rhs.m_matrixData[i];
		
	iMatrix<T> result(numRows, numCols, tempResult);
	delete[] tempResult;
	return result;
}


template <class T>
iMatrix<T> operator* (const iMatrix<T>& lhs, const T& rhs)
{
	int numRows = lhs.m_nRows;
	int numCols = lhs.m_nCols;
	int numElements = lhs.m_matrixData.size();
	T *tempResult = new T[numElements];
	for (int i=0; i<numElements; ++i)
		tempResult[i] = lhs.m_matrixData[i] * rhs;
		
	iMatrix<T> result(numRows, numCols, tempResult);
	delete[] tempResult;
	return result;
}

*/
// matrix * scaler
template <class T>
iMatrix<T> operator* (const iMatrix<T>& lhs, const T& rhs)
{
    iMatrix<T> result(lhs.m_nRows, lhs.m_nCols);

    std::transform(
        lhs.m_matrixData.begin(), 
        lhs.m_matrixData.end(),   
        result.m_matrixData.begin(), 
        [&rhs](const T& element) { return element * rhs; }
    );

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
		iMatrix<T> rhs_T(r_numCols, r_numRows);
		for (int i = 0; i < r_numRows; ++i) {
			for (int j = 0; j < r_numCols; ++j) {
				rhs_T(j, i) = rhs(i, j);
			}
		}

		T *tempResult = new T[lhs.m_nRows * rhs.m_nCols];

		for (int lhsRow=0; lhsRow<l_numRows; lhsRow++) 
		{
			for (int rhsCol=0; rhsCol<r_numCols; rhsCol++) 
			{
				T elementResult = static_cast<T>(0.0);
				int lhsLinearIndex = (lhsRow * l_numCols); 
				int rhsTLinearIndex = (rhsCol * r_numRows); 
				
				
				for (int k=0; k<l_numCols; k++) 
				{
					
					elementResult += (lhs.m_matrixData[lhsLinearIndex + k] * rhs_T.m_matrixData[rhsTLinearIndex + k]);
				}
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
	
	if (this != &rhs)
	{
		m_nRows = rhs.m_nRows;
		m_nCols = rhs.m_nCols;
		m_matrixData = rhs.m_matrixData; 
	}
	
	return *this;
}


template <class T>
iMatrix<T> operator* (const std::vector<T>& v1, const std::vector<T>& v2)
{
    int N = v1.size();
    int M = v2.size(); 
    
    iMatrix<T> result(N, M); 
    for (int i = 0; i < N; ++i) 
    {
        T v1_i = v1[i]; 
        
        for (int j = 0; j < M; ++j) 
        {
            result(i, j) = v1_i * v2[j]; 
        }
    }

    
    return result;
}

template <class T>
iMatrix<T>& iMatrix<T>::operator+= (const iMatrix<T>& rhs)
{

    if (m_nRows != rhs.m_nRows || m_nCols != rhs.m_nCols)
    {
        throw std::invalid_argument("Error: Las matrices deben tener las mismas dimensiones para la suma in-place.");
    }

    int numElements = m_matrixData.size();
    T* data_lhs = m_matrixData.data();
    const T* data_rhs = rhs.m_matrixData.data(); 
    for (int i = 0; i < numElements; ++i)
    {
        data_lhs[i] += data_rhs[i]; 
    }
    return *this;
}

template <class T>
iMatrix<T>& iMatrix<T>::operator*= (const T& rhs)
{
    int numElements = m_matrixData.size();
    T* data_lhs = m_matrixData.data();
    for (int i = 0; i < numElements; ++i)
    {
        data_lhs[i] *= rhs; 
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

// Function to test whether the matrix is square.
template <class T>
bool iMatrix<T>::IsSquare()
{
	if (m_nCols == m_nRows)
		return true;
	else
		return false;
}

template <class T>
T iMatrix<T>::Determinant() const
{
    if (!this->IsSquare()){
        return static_cast<T>(0.0);
    }

    int N = this->m_nRows;
    if (N == 1){
        return (*this)(0, 0); 
    }
    iMatrix<T> tempMatrix = *this; 

    T det = static_cast<T>(1.0);
    int sign = 1;

    for (int k = 0; k < N; ++k){
        int pivotRow = k;
        T maxValue = std::abs(tempMatrix(k, k));
        
        for (int i = k + 1; i < N; ++i) 
        {
            if (std::abs(tempMatrix(i, k)) > maxValue) 
            {
                maxValue = std::abs(tempMatrix(i, k));
                pivotRow = i;
            }
        }

        if (pivotRow != k){
            for (int j = 0; j < N; ++j){
                std::swap(tempMatrix(k, j), tempMatrix(pivotRow, j)); 
            }

            sign = -sign; 
        }

        if (std::abs(tempMatrix(k, k)) < static_cast<T>(1e-9)){
            return static_cast<T>(0.0); 
        }

        T pivotValue = tempMatrix(k, k);
        for (int i = k + 1; i < N; ++i) {
            T factor = tempMatrix(i, k) / pivotValue;
            
            for (int j = k; j < N; ++j){
                tempMatrix(i, j) -= factor * tempMatrix(k, j);
            }
        }

        det *= pivotValue;
    }
    return sign * det;
}

template <class T>
std::vector<T> iMatrix<T>::Vectorize() const
{
    return m_matrixData;
}

template <class T>
iMatrix<T> iMatrix<T>::Transpose() const
{
    iMatrix<T> result(m_nCols, m_nRows);  // ← filas y columnas intercambiadas

    // Recorremos todos los elementos de la matriz original
    for (int i = 0; i < m_nRows; ++i)
    {
        for (int j = 0; j < m_nCols; ++j)
        {
            // Elemento (i,j) en la original → posición (j,i) en la transpuesta
            result(j, i) = (*this)(i, j);
        }
    }

    return result;
}



template <class T>
T dot_product (const std::vector<T>& v1, const std::vector<T>& v2)
{
    if (v1.size() != v2.size())
    {
        throw std::invalid_argument("Error: Los vectores deben tener la misma dimensión para el Producto Punto.");
    }
    
    int N = v1.size();
    T resultValue = static_cast<T>(0.0);
    
    for (int i = 0; i < N; ++i)
    {
        resultValue += v1[i] * v2[i];
    }

    return resultValue;
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



#endif