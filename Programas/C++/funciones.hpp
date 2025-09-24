/******************************************************************************* 
 *AUTHOR: Carlos Palomera Oliva
 *DATE: MARCH 2025
 *******************************************************************************/
#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>  
#include <fstream>
#include <string>
#include <limits>
#include "iMatrix.cpp"
#include "C:\Librerias\eigen-3.4.0\Eigen\Sparse"
#include "C:\Librerias\eigen-3.4.0\Eigen\Dense"    
#ifndef FUNCIONES_HPP
#define FUNCIONES_HPP




/**************************************************************************************
 * Structure of the Data File for Matrix
 * Nrow Ncol
 * vector of size Ncol representing the first variable where elements are separated by space Ex(x1 x2 x3 ...)
 * vector of size Ncol ...
 * .
 * .  Nrow times in total
 * . 
 * 
 **************************************************************************************/
/*
*Input: Direction of the data file
*Output: iMatrix filled with the control points in data file
*/
iMatrix<double> ReadDataFile_Matrix(std::string direction);

/**************************************************************************************
 * Structure of the Data File for Tensors
 * Ncol NSplits   *
 * vector of size Ncol representing the first variable where elements are separated by space Ex(x1 x2 x3 ...)
 * vector of size Ncol ...
 * .
 * .  3 times in total
 * . 
 * --------
 * The previous repeated NSplits totals
 **************************************************************************************/
/*
*Input: Direction of the data file
*Output: vector of iMatrixes filled with the control points in data file
*/
std::vector<iMatrix<double>> ReadDataFile_CtrlPts(std::string direction);


/*
*Input: Value of the point where the polinomial is evaluated t, and degree of the polinomial n
*Condition: t€[0,1]
*Output: Values of all the Bernstein polinomial of degree n evaluated at t
*/
std::vector<double> BernsteinPolinomial(double t, int n);


/*
*Input: iMatrix of control points and parameter of the curve t
*Conditions: t€[0,1]
*Output: Value of the Bezier curve of such control points in the point t
*/
std::vector<double> BezierCurve_Point_t(const iMatrix<double> &PtosControl, double t);


/*
*Input: iMatrix of control points and number of steps for the curve
*Conditions: Steps must be >= 1
*Output: Value of the Bezier curve at every point from 0 to 1 separated in N_Steps steps
*/
iMatrix<double> BezierCurve(const iMatrix<double> &PtosControl, double N_Steps);

/*
*Input: Tensor of control points and the 2 parameters of the curve t1, t2
*Conditions: t1,t2€[0,1]
*Output: Value of the Bezier surface at of such control points in the points t1 t2
*/
std::vector<double> BezierSurface_Point_t1_t2(const std::vector<iMatrix<double>> &PtosControl, double t1, double t2);

/*
*Input: iMatrix of control points and number of steps for the curve in both parameters
*Conditions: N_Steps must be >= 1
*Output: Value of the Bezier surface at every point from 0 to 1 separated in N_Steps along both parameters
*/
std::vector<iMatrix<double>> BezierSurface(const std::vector<iMatrix<double>> &PtosControl, double N_Steps);

/*
*Input: Bezier surface in form of a tensor
*Output: Writes a file of 3 matrices one for each coordinate
*/
void Write_BezierSurface(const std::vector<iMatrix<double>> &Tensor, std::string name);

/*
*Input: CtrlPts in format described before
*Output: Writes a file of 3 matrices one for each coordinate
*/
void Write_CtrlPts(const std::vector<iMatrix<double>> &CtrlPts, std::string name);



/************************************************************************************
*
* Classic funtions from the Nurbs Book
*
************************************************************************************/

/*
*Input: A KnotVector,t the point of evalution, its polynomial degree, and n=the size of the Nurbs vector-p-2
*Output: the span index
*/
int FindSpan(const std::vector<double> &KnotVector, double t, int p, int n);

/*
*Input: i span index, t eval point, p polinomial degree and KnotVector
*Output: vector of basis funcs
*/
std::vector<double> BasisFuns(int i, double t, int p,const std::vector<double> &KnotVector);

/*
*Input: i span index, t eval point, p polinomial degree and KnotVector
*Output: iMatrix where the i-th row is the i-th degree basis fun
*/
iMatrix<double> AllBasisFuns(int i, double t, int p,const std::vector<double> &KnotVector);

/*
*Input: i span index, t eval point, p polinomial degree, k derivate degree, KnotVector and an output vector ders;
*Output: value of the k-derivate at point t in the ders vector on the k-row of the matrix
*/
iMatrix<double> DerBasisFuns(int span, double t, int p, int k,const std::vector<double> &KnotVector);

/*
*Input: matrix of BasisFunctions or its derivates and the path where they will be written
*Output: Writes the Functions as a matrix on selected path
*/
void Write_BasisFunctions(const iMatrix<double> &BasisFunctions, std::string name);

/************************************************************************************
*
* Classic funtions for p spline curves/surfaces
*
************************************************************************************/

/*
*Input: n=the size of the Nurbs vector-p-2, p polynomial degree, KnotVector, the CtrlPts ,t the point of evalution
*Output: Evaluation of the curve at point t
*/
std::vector<double> CurvePoint(int n, int p,const std::vector<double> &KnotVector,const iMatrix<double> &CtrlPts, double t);

/*
*Input: n=the size of the Nurbs vector-p-2, p polynomial degree, KnotVector, the CtrlPts ,t the point of evalution and d highest derivate order to compute
*Output: iMatrix where the i-th row is the evaluation of the i-th derivate
*/
iMatrix<double> CurveDerivsAlg1(int n, int p,const std::vector<double> &KnotVector,const iMatrix<double> &CtrlPts, double t, int d);

/*
*Input: n=the size of the Nurbs vector-p-2, p polynomial degree, KnotVector, the CtrlPts ,t the point of evalution, d highest derivate order to compute, r1 and r2 specify the intervals of CtrlPoints
*Output: iMatrix where the i-th row is the vector of ctrl points
*/
iMatrix<std::vector<double>> CurveDerivsCpts(int n, int p,const std::vector<double> &KnotVector,const iMatrix<double> &CtrlPts,  int d, int r1, int r2);

/*
*Input: n=the size of the Nurbs vector-p-2, p polynomial degree, KnotVector, the CtrlPts ,t the point of evalution and d highest derivate order to compute
*Output: iMatrix where the i-th row is the evaluation of the i-th derivate
*/
iMatrix<double> CurveDerivsAlg2(int n, int p,const std::vector<double> &KnotVector,const iMatrix<double> &CtrlPts, double t, int d);

/*
*Input: Curve and path
*Output: Writes the values of the Curve in the selected path
*/
void Write_Curve(const iMatrix<double> &Curve, std::string filename);

/*
*Input: n_i=the size of the Nurbs vector_i-p-2, p_i polynomial degree, KnotVector_i, the CtrlPts ,t_i the point of evalution
*Output: Evaluation of the surface at point (t1,t2)
*/
std::vector<double> SurfacePoint(int n1, int p1,const std::vector<double> &KnotVector1, int n2, int p2,const std::vector<double> &KnotVector2,const std::vector<iMatrix<double>> &CtrlPts, double t1, double t2);

/*
*Input: n_i=the size of the Nurbs vector_i-p-2, p_i polynomial degree, KnotVector_i, the CtrlPts ,t_i the point of evalution and d maximum order of the derivates
*Output: Evaluation of the derivates up to degree d (maximum is p or q and would overwrite it if d>p,q) in (t1,t2) where Output[k_1][k_2] is been derived k_i times respect to the ith variable 
*/
iMatrix<std::vector<double>> SurfaceDerivsAlg1(int n1, int p1,const std::vector<double> &KnotVector1, int n2, int p2,const std::vector<double> &KnotVector2,const std::vector<iMatrix<double>> &CtrlPts, double t1, double t2, int d);


/************************************************************************************
*
* Classic and support funtions for NURBS curves/surfaces
*
************************************************************************************/


/*
*Input: Vector of CtrlPts and their weights at each point in a separated vector
*Output: iMatrix that each column is each of the coordinates of each CtrlPt except its last component that is it's respective weight
*/
iMatrix<double> WeightCtrlPts(const iMatrix<double> &CtrlPts,const std::vector<double> &Weights);

/*
*Input: CtrlPts multiplied by their Weight with the respective weight at the end of each column
*Output: The original CtrlPts
*/
iMatrix<double> UnWeightCtrlPts(const iMatrix<double> &CtrlPtsW);

/*
*Input: Vector of CtrlPts and their weights at each point in a separated vector
*Output: std::vector<iMatrix<double>> that keeps the format of the CtrlPts for Surface operations with and extra 4th dimension on each coord thats's occupied by the weight of such pt 
*/
std::vector<iMatrix<double>> WeightCtrlPtsSurf(const std::vector<iMatrix<double>> &CtrlPts,const iMatrix<double> &Weights);

/*
*Input: CtrlPts multiplied by their Weight with the respective weight at the end of each column
*Output: The original CtrlPts
*/
std::vector<iMatrix<double>> UnWeightCtrlPtsSurf(const std::vector<iMatrix<double>> &CtrlPtsW);

/*
*Input: n=the size of the Nurbs vector-p-2, p polynomial degree, KnotVector, the CtrlPts weighted with theris weight in the last component ,t the point of evalution
*Output: Evaluation of the NURBS curve at point t
*/
std::vector<double> CurvePointRational(int n, int p,const std::vector<double> &KnotVector,const iMatrix<double> &CtrlPtsWeighted, double t);

/*
*Input: binomial coeficients k, i
*Output: combinatiorial number given by upper number k and lower number i
*/
int Binomial(int k, int i);

/*
*Input: i span index, t eval point, p polinomial degree and KnotVector, weights vector
*Output: vector of rational basis funcs
*/
std::vector<double> RationalBasisFuns(int i, double t, int p,const std::vector<double> &KnotVector,const std::vector<double> &WeightVector);

/*
*Input: i span index, t eval point, k derivate degree, p polinomial degree and KnotVector, weights vector
*Output: vector of rational basis funcs
*/
iMatrix<double> RatDersBasisFuns(int i, double t, int p, int k, const std::vector<double> &KnotVector,const std::vector<double> &WeightVector);

/*
*Input: Aders: derivates of the CtrlPtsWeighted using CurvDersAlg1 evaluated at some t, wders: same as Aders but with the weighted vector, maximum degree of the derivate
*              -Aders = CurveDerivsAlgi(n, p,KnotVector, WeightedCtrlPts, t, d);
*              -wders = CurveDerivsAlgi(n, p,KnotVector, WeightVector, t, d);  
*Output: iMatrix on which each row is the i-th derivate of the rational curve made of the weighted CtrlPts and the weight vector
*/
iMatrix<double> RatCurveDerivs(const iMatrix<double> &Aders,const iMatrix<double> &wders, int d);

/*
*Input: n_i=the size of the Nurbs vector_i-p-2, p_i polynomial degree, KnotVector_i, the CtrlPtsWeighted ,t_i the point of evalution
*Output: Evaluation of the surface at point (t1,t2)
*/
std::vector<double> SurfacePointRational(int n1, int p1,const std::vector<double> &KnotVector1, int n2, int p2,const std::vector<double> &KnotVector2,const std::vector<iMatrix<double>> &CtrlPtsWeighted, double t1, double t2);

/*
*Input: derivates of the CtrlPtsWeighted using CurvDersAlg1 evaluated at some t, wders: same as Aders but with the weighted vector, maximum degree of the derivate
*Output: iMatrix on which each row is the i-th derivate of the rational surface made of the weighted CtrlPts and the weight vector
*/
iMatrix<std::vector<double>> RatSurfaceDerivs(const iMatrix<std::vector<double>> &Aders,const iMatrix<double> &wders,int d);

/************************************************************************************
*
* Fundamental Geometric Algorithms
*
************************************************************************************/

/*
*Input: KnotVector and element you want to introduce
*Output: index i such that KnotVector[i-1] <= e < KnotVector[i]
*Note: Applying insertion with this index using KnotVector.insert(KnotVector.begin()+index+1, r, e) will yield the KnotVector we want (r is number of insertions)
*/
int findInsertionIndex(const std::vector<double> &KnotVector, double e);

/*
*Input: ordered vector and element you want to check it's number of ocurrences
*Output: number of times e appears on the vector
*/
int countOccurrences(const std::vector<double>& v, double e);

/*
*Input: A KnotVector,t the point of evalution, its polynomial degree, and n=the size of the Nurbs vector-p-2
*Output: the span index at k, its multiplicity at s
*/
void FindSpanMult(const std::vector<double> &KnotVector, double t, int p, int n,int &k, int &s);

/*
*Input: np number of initial CtrlPts-1, p polinomial degree, the initial KnotVector, the initial CtrlPts Weighted, u value of the inserted knot,
        index such that KnotVector[k]<= u < KnotVector[k+1], s initial multiplicity of u in KnotVector, r number of insertions of u in Knotvector 
*Output: New CtrlPts Weighted that presserves the initial curve
*Note: !! KnotVector will be updated !!
*/
iMatrix<double> CurveKnotsIns(int np, int p, std::vector<double> &KnotVector,const iMatrix<double> &CtrlPtsW, double u, int index, int s, int r);

/*
*Input: n=the size of the Knot vector-p-2, p polynomial degree, KnotVector, CtrlPts Weigthed, u point of evaluation 
*Output: Evaluation of the curve given by the KnotVector and the ControlPts at the given point of evaluation u
*/
std::vector<double> CurvePntByCornerCut(int n, int p, std::vector<double> &KnotVector, iMatrix<double> &CtrlPtsW, double u);

/*
*Input: npi size of KnotVectori-pi-2, pi degree of KnotVectori, KnotVectori for i=1,2, CtrlPts Weigthed, dir(true= knot inserted in Knotvector1, false= knot inserted in Knotvector2
        uv knot we want to insert, index such that KnotVectori[k]<= u < KnotVectori[k+1], s initial multiplicity of u in KnotVectori, r number of insertions of u in Knotvectori
*Output: the KnotVectori will be updated and it will return new CtrlPtsW such that they generate the same surface.
*/
std::vector<iMatrix<double>> SurfaceKnotIns(int np1, int p1, std::vector<double> &KnotVector1, int np2, int p2, std::vector<double> &KnotVector2,const std::vector<iMatrix<double>> &CtrlPtsW, bool dir,
                    double uv, int index, int r, int s);

/*
*Input: n=Knotvector.size()-2-p, p polynomial degree, KnotVector, CtrlPtsW, Insertions: vector of new elements to insert into KnotVector, r=Insertion.size()-1. 
*Output: Update of the Knotvector with the new inserions, CtrlPtsW that generate the same curve with the new KnotVector.
*/
iMatrix<double> RefineKnotVectCurve(int n, int p, std::vector<double> &KnotVector,const iMatrix<double> &CtrlPtsW,const std::vector<double> &Insertions, int r);

/*
*Input: npi size of KnotVectori-1, pi degree of KnotVectori, KnotVectori for i=1,2, CtrlPts Weigthed, dir(true= knot inserted in Knotvector1, false= knot inserted in Knotvector2
        Insertions: vector of new elements to insert into KnotVector,  r=Insertion.size()-1.
*Output: the KnotVectori will be updated and it will return new CtrlPtsW such that they generate the same surface.
*/
std::vector<iMatrix<double>> RefineKnotVectSurface(int np1, int p1, std::vector<double> &KnotVector1, int np2, int p2, std::vector<double> &KnotVector2,const std::vector<iMatrix<double>> &CtrlPtsW, bool dir,
                    const std::vector<double> Insertions, int r);

/************************************************************************************
*
* Functions for the Assembly of the Matrix
*
************************************************************************************/

/*
*Input: n: nº pts of integration
*Output: nodes will contain the roots of the legendre polynomial and weights its weights (all in double precision).
*/
void legendre_pol(int n, Eigen::VectorXd &nodes, Eigen::VectorXd &weights);

/*
*Input: ordered vector
*Output: the original vector without duplicates
*/
std::vector<double> removeDuplicates(const std::vector<double>& input);

/*
*Input: number of elements of evaluation
*Output: vector of nElements+1 uniformly distributed
*/
std::vector<double> createSubIntervals(const double nElements, const double lower_limit, const double upper_limit);

/*
*Input: v1 vector where we remove the elements that appear in v2, v2 KnotVector
*Output: v1 without the elements in v2
*/
std::vector<double> removeMatches(const std::vector<double>& v1,  const std::vector<double>& v2);

/*
*Input: n=the size of the Knot vector-p-2, i span index, t eval point, p polinomial degree, number of points of evaluations in the integral,
*        KnotVector, upper and lower limits of the integration, weights vector, nodes of the legendre polynomial,  CtrlPtsW the control points weighted, and the function f.
*Output: BasisFunsEvals and DerBasisFunsEvals each row will store the evaluation of the value or derivate of the basis fun alongside each evaluation point, same for the JacobianEvals, InverseJacobianEvals and funcEvals.
*/
void D1_element_eval(int n, int i, int p, int nEvals, double lower_limit, double upper_limit,const std::vector<double> &KnotVector,const std::vector<double> &Weights,
                     const Eigen::VectorXd &nodes,const iMatrix<double> &CtrlPtsW, std::function<double(double)> f,
                     iMatrix<double> &BasisFunsEvals, iMatrix<double> &DerBasisFunsEvals,std::vector<double>& JacobianEvals, std::vector<double>& InverseJacobianEvals, std::vector<double>& funcEvals, double Lenght=1.0);

/*
*Input: p polinomial degree, upper and lower limits of the integration, number of points of evaluations in the integral, evaluations of RationalDerivates, evaluation of the Inverse Jacobian and the legendre weights.
*Output: Matrix of all the permutations of evaluations where the elements i,j of the matrix correspond to the BasisFuns(init+i)*BasisFuns(init+j) where init=spanIndex()
*/                        
iMatrix<double> gauss_legendre_cuadrature_integral_bilinealForm(int p, double lower_limit, double upper_limit, int nEvals,const iMatrix<double> &DerBasisFunsEvals,const std::vector<double>& InverseJacobianEvals,const Eigen::VectorXd &weights);

/*
*Input: p polinomial degree, upper and lower limits of the integration, number of points of evaluations in the integral, evaluations of RationalBasisFuns, evaluations of the Jacobian, evaluations of the function and the legendre weights.
*Output: Vector of evaluations which its components are l(index) where init=spanIndex()
*/                        
std::vector<double> gauss_legendre_cuadrature_integral_linealForm(int p, double lower_limit, double upper_limit, int nEvals,const iMatrix<double> &BasisFunsEvals,const std::vector<double>& JacobianEvals,const std::vector<double>& funcEvals,const Eigen::VectorXd &weights);

/*
*Input: p polinomial degree, size=LineraForm.size()-1, SparseMatrix with the system already assambled, same with the LinearForm, Modify0(1) true if we have dirichlet condition at the start(end) of the boundary,
*       vValueAt0(1) value set for the boundary condition, default is set at 0.
*Output: Modification of the sparse matrix and the linear form to keep the system as symetric as posible.
*/  
void impose_Dirichlet_Condition(int p, int size, Eigen::SparseMatrix<double> &global, Eigen::VectorXd &LinearForm, bool Modify0, bool Modify1, double ValueAt0=0.0, double ValueAt1=0.0);

/*
*Input: size=LineraForm.size()-1, LinearForm, Modify0(1) true if we have dirichlet condition at the start(end) of the boundary,
*       vValueAt0(1) value set for the boundary condition, default is set at 0.
*Output: Modification of the linear form applying the conditions to the borders.
*/ 
void impose_Newmann_Condition(int size, Eigen::VectorXd &LinearForm, bool Modify0, bool Modify1, double ValueAt0=0.0, double ValueAt1=0.0);

/*
*Input: p polinomial degree, size=LineraForm.size()-1, SparseMatrix with the system already assambled, same with the LinearForm, Modify0(1) true if we have dirichlet condition at the start(end) of the boundary,
*       vValueAt0(1) value set for the boundary condition, default is set at 0.
*       values alpha and beta come from the expresion alpha*u(x)+beta*u'(x)=ValueAtx.
*Output: Modification of the sparse matrix and the linear form to keep the system as symetric as posible.
*/  
void impose_Robin_Condition(int p, int size, Eigen::SparseMatrix<double> &global, Eigen::VectorXd &LinearForm, double alpha, double beta, bool Modify0, bool Modify1, double ValueAt0=0.0, double ValueAt1=0.0);

/*
*Input: vector of evaluations of the function, name of hte file where we write the evaluations
*Output: writes the evaluation of the function in a txt
*/  
void writeNumericFunction(Eigen::VectorXd funEvals,const std::vector<double> &KnotVector,const std::vector<double> &WeightVector, int n, int p, int nEvals, double lowerLimit, double upperLimit,  std::string name,
                            Eigen::VectorXd nodes, Eigen::VectorXd weights, std::vector<double> time, int nElements); 

/*
*Input: the function we evaluate, lower and upper limits of evaluation and number of steps of the evaluation, name of hte file where we write the evaluations,
        for the evaluation of the curve point, CtrlPtsW, KnotVector, n = KnotVector.size()-p-2, p polynomial degree.
*Output: writes the evaluation of the function in a txt
*/  
void writeAnalyticFunction(std::function<double(double)> f,std::function<double(double)> f_Ders, double lowerLimit, double upperLimit, double nEvals,const iMatrix<double> &CtrlPtsW,const std::vector<double> &KnotVector,
                        int n, int p, std::string name,Eigen::VectorXd nodes, Eigen::VectorXd weights, std::vector<double> time, double Lenght);

/*
*Input: vector of evaluations of the function, name of hte file where we write the evaluations
*Output: writes the evaluation of the derivates of the function in a txt
*/ 
/*
void writeFunctionDers(Eigen::VectorXd funEvals,const std::vector<double> &KnotVector,const iMatrix<double> &CtrlPtsW,const std::vector<double> &WeightVector, int n, int p, int nEvals, double lowerLimit, double upperLimit,  std::string name);
*/
/*
*Input: the function we evaluate, lower and upper limits of evaluation and number of steps of the evaluation, name of hte file where we write the evaluations,
        for the evaluation of the curve point, CtrlPtsW, KnotVector, n = KnotVector.size()-p-2, p polynomial degree.
*Output: writes the evaluation of the function in a txt
*/ 
/* 
void writeFunctionDers(std::function<double(double)> f, double lowerLimit, double upperLimit, double nEvals,const iMatrix<double> &CtrlPtsW,const std::vector<double> &KnotVector,
                    int n, int p, std::string name, double Lenght=1.0);
*/
/*
*Input: lowerLimit of integration, upperLimit of integration, number of subintervals for the quadrature, CtrlPtsW, Knotvector, n, p
*Output: numerical value of the integral between lowerLimit and upperLimit of the norm of the derivate of the curve
*/ 
double IntegrateNormDer(double lowerLimit, double upperLimit, int nEvals,const iMatrix<double> &CtrlPtsW,const std::vector<double> &KnotVector, int n, int p);

/*
*Input: e fisical parameter, lowerLimit of integration, upperLimit of integration, number of subintervals for the quadrature, CtrlPtsW, Knotvector, n, p
*Output: numerical value parametric space parameter such that the analytic solution of it its equal to the fisical solution at e
*/ 
inline double Fisical_to_Parametric(double e, double lowerLimit, double Lenght, int nEvals,const iMatrix<double> &CtrlPtsW,const std::vector<double> &KnotVector, int n, int p){
        return IntegrateNormDer(lowerLimit,e,nEvals, CtrlPtsW, KnotVector, n, p)/Lenght;
}

#endif
