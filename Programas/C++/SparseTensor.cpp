#include <vector>
#include "C:\Librerias\eigen-3.4.0\Eigen\Sparse"
#include <unordered_map>
#include <functional>
#include <stdexcept>
#include <algorithm>

struct VectorHash {
    std::size_t operator()(const std::vector<int>& v) const {
        std::size_t seed = 0;
        for (int i : v) {
            seed ^= std::hash<int>{}(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

template <class T>
class SparseTensor {
private:
    std::unordered_map<std::vector<int>, T, VectorHash> data;
    int dim;

public:
    SparseTensor(int dimension) : dim(dimension) {}

    // Fijar un valor
    void SetElement(const std::vector<int>& index, T value) {
        if (index.size() != dim) {
            throw std::invalid_argument("Index dimension mismatch.");
        }
        if (value != 0.0) {
            data[index] = value;
        } else {
            data.erase(index); // eliminar si es 0
        }
    }

    // Leer un valor
    T GetElement(const std::vector<int>& index) const {
        if (index.size() != dim) {
            throw std::invalid_argument("Index dimension mismatch.");
        }
        auto it = data.find(index);
        if (it != data.end()) {
            return it->second;
        }
        if constexpr (std::is_arithmetic_v<T>) {
		// Si T es un tipo numérico (int, float, double, etc.)
		    return  0.0f;
        } 
        else {
            // Si T es un contenedor (std::vector<float>, etc.)
            throw std::invalid_argument("No default value.");
        }
    }

    // Verificar si existe un valor no cero
    bool exists(const std::vector<int>& index) const {
        return data.find(index) != data.end();
    }

    // Iterar sobre elementos no nulos
    void forEach(std::function<void(const std::vector<int>&, T)> func) const {
        for (const auto& pair : data) {
            func(pair.first, pair.second);
        }
    }
 
    void forEachComp(int comp, int targetCoord, std::function<void(const std::vector<int>&, double)> func) const {
        if (comp < 0 || comp >= dim) {
            throw std::out_of_range("Component index out of range");
        }
        for (const auto& [idx, val] : data) {
            if (comp >= 0 && comp < idx.size() && idx[comp] == targetCoord) {
                func(idx, val);
            }
        }
    }

    int dimension() const {
        return dim;
    }

    /************************************************************************************
    *
    * Classic and support funtions for NURBS curves/surfaces
    *
    ************************************************************************************/

    int flattenIndex(const std::vector<int>& idx, const std::vector<int>& dims) const {
        int flat = 0;
        int stride = 1;
        for (int i = idx.size() - 1; i >= 0; --i) {
            flat += idx[i] * stride;
            stride *= dims[i];
        }
        return flat;
    }

    using Triplet = Eigen::Triplet<T>;


    Eigen::SparseMatrix<T> toEigenSparse() const {
        if (dimension() < 2){
            throw std::invalid_argument("Tensor dimension must be >= 2");
        }
        
        std::vector<int> maxDims(dim, 0);

        // Encontrar tamaño máximo por dimensión
        forEach([&](const std::vector<int>& idx, T val) {
            for (int i = 0; i < dim; ++i)
                if (idx[i] > maxDims[i]) maxDims[i] = idx[i];
        });

        for (auto& d : maxDims) ++d; // Añadir 1 para cubrir el índice

        // Producto cartesiano para filas y columnas
        int rows = 1, cols = 1;
        for (int i = 0; i < dim / 2; ++i) rows *= maxDims[i];
        for (int i = dim / 2; i < dim; ++i) cols *= maxDims[i];

        std::vector<Triplet> triplets;

        forEach([&](const std::vector<int>& idx, T val) {
            int row = flattenIndex(
                std::vector<int>(idx.begin(), idx.begin() + dim / 2),
                std::vector<int>(maxDims.begin(), maxDims.begin() + dim / 2));
            int col = flattenIndex(
                std::vector<int>(idx.begin() + dim / 2, idx.end()),
                std::vector<int>(maxDims.begin() + dim / 2, maxDims.end()));

            triplets.emplace_back(row, col, val);
        });

        Eigen::SparseMatrix<T> mat(rows, cols);
        mat.setFromTriplets(triplets.begin(), triplets.end());
        return mat;
    }   

};