bool PrintMatrixes=false;
bool PrintMatrixes2=false;
bool ShowValues=false;
bool ShowIterations=false;
bool ShowIterations2=false;
bool ShowConditions=false;
bool ShowSolutions=false;

//Initial definitions of variables
int pY= 2;
int pX= 2;
int hY= 512;
int hX= 512;
int nEvalsY=pY+1;
int nEvalsX=pX+1;
int subMatSize=nEvalsX*nEvalsY;
int nElementsY=hY+pY;
int nElementsX=hX+pX;
int size = nElementsY*nElementsX;

Eigen::SparseMatrix<double> global(size, size);
Eigen::VectorXd LinearForm(size);
std::vector<double> KnotVectorY = {0,0,0,1,1,1};
std::vector<double> KnotVectorX = {0,0,0,1,1,1};
        
double lower_limitY=KnotVectorY[0];
double upper_limitY=KnotVectorY[KnotVectorY.size()-1];
double lower_limitX=KnotVectorX[0];
double upper_limitX=KnotVectorX[KnotVectorX.size()-1];



auto analyticSol = [](double x, double y){
    //double pi=3.141592653589793;
    return (1-x*x-y*y)/4;
};

auto analyticSol_dx = [](double x, double y) {
    //double pi = 3.141592653589793;
    return -2*x/4;
};

auto analyticSol_dy = [](double x, double y) {
    //double pi = 3.141592653589793;
    return -2*y/4;
};

auto f = [](double x, double y){
    //double pi=3.141592653589793;
    return 1;
};

std::vector<double> data={1,std::sqrt(2)/2,1,
                            std::sqrt(2)/2,1,std::sqrt(2)/2,
                            1,std::sqrt(2)/2,1,
                            
                            };

iMatrix<double> Weights(3,3, data.data());
std::vector<iMatrix<double>> CtrlPts = ReadDataFile_CtrlPts(PATH + "EjSuperficieCircle.txt");
std::vector<iMatrix<double>> CtrlPtsW=WeightCtrlPtsSurf(CtrlPts,Weights);
Eigen::VectorXd nodesY, nodesX, weightsY, weightsX;
legendre_pol(nEvalsY, nodesY, weightsY);
legendre_pol(nEvalsX, nodesX, weightsX);


if(ShowValues){
    std::cout<< "---------------KnotVectorY-----------------" << std::endl;
for(unsigned int i=0; i<KnotVectorY.size(); i++){
    std::cout << KnotVectorY[i] << ", ";
}
std::cout << std::endl;

std::cout<< "---------------KnotVectorX-----------------" << std::endl;
for(unsigned int i=0; i<KnotVectorX.size(); i++){
    std::cout << KnotVectorX[i] << ", ";
}
std::cout << std::endl;


std::cout << "lower_limitY= " << lower_limitY << std::endl;
std::cout << "upper_limitY= " << upper_limitY << std::endl;
std::cout << "lower_limit2= " << lower_limitX << std::endl;
std::cout << "upper_limit2= " << upper_limitX << std::endl;

std::cout<<"----CtrlPts----" << std::endl;
for(unsigned int i=0; i<CtrlPts.size();i++){
    CtrlPts[i].PrintMatrix();
    std::cout<<"-------------------" << std::endl;
}
std::cout<<"---Matriz_Pesos---" << std::endl;
Weights.PrintMatrix();
std::cout<<"-----CtrlPtsW-----" << std::endl;

for(unsigned int i=0; i<CtrlPtsW.size();i++){
    CtrlPtsW[i].PrintMatrix();
    std::cout<<"-------------------" << std::endl;
}
}

std::vector<double> InsertionsY=createSubIntervals(hY, lower_limitY, upper_limitY);
InsertionsY=removeMatches(InsertionsY, KnotVectorY);


std::vector<double> InsertionsX=createSubIntervals(hX, lower_limitX, upper_limitX);
InsertionsX=removeMatches(InsertionsX, KnotVectorX);



std::vector<iMatrix<double>> NewCtrlPtsWintermedieate= RefineKnotVectSurface(KnotVectorY.size()-2-pY,pY,KnotVectorY,KnotVectorX.size()-2-pX,pX,
                                                                            KnotVectorX, CtrlPtsW, true, InsertionsY, InsertionsY.size()-1);
std::vector<iMatrix<double>> NewCtrlPtsW= RefineKnotVectSurface(KnotVectorY.size()-2-pY,pY,KnotVectorY,KnotVectorX.size()-2-pX,pX,
                                                                KnotVectorX, NewCtrlPtsWintermedieate, false, InsertionsX, InsertionsX.size()-1);


NewCtrlPtsWintermedieate.clear();

iMatrix<double> NewWeights(NewCtrlPtsW.size(),NewCtrlPtsW[0].GetNumCols());
for(unsigned int i=0; i<NewWeights.GetNumRows(); i++){
    for(unsigned int j=0; j<NewWeights.GetNumCols(); j++){
        NewWeights(i,j)=NewCtrlPtsW[i](3,j);
    }
}
std::vector<double> IntervalsY=removeDuplicates(KnotVectorY);
std::vector<double> IntervalsX=removeDuplicates(KnotVectorX);





if(ShowValues){

    std::cout<< "---------------InsertionsY-----------------" << std::endl;
    for(unsigned int i=0; i<InsertionsY.size(); i++){
        std::cout << InsertionsY[i] << ", ";
    }
    std::cout << std::endl;

    std::cout<< "---------------InsertionsX-----------------" << std::endl;
    for(unsigned int i=0; i<InsertionsX.size(); i++){
        std::cout << InsertionsX[i] << ", ";
    }
    std::cout << std::endl;


    std::cout<<"------NewCtrlPtsW-----" << std::endl;
    for(unsigned int i=0; i<NewCtrlPtsW.size();i++){
        NewCtrlPtsW[i].PrintMatrix();
        std::cout<<"-------------------" << std::endl;
    }

    std::cout<<"------NewWeights-----" << std::endl;
    NewWeights.PrintMatrix();

    std::cout<< "---------------NewKnotVectorY-----------------" << std::endl;
    for(unsigned int i=0; i<KnotVectorY.size(); i++){
        std::cout << KnotVectorY[i] << ", ";
    }
    std::cout << std::endl;

    std::cout<< "---------------NewKnotVectorX-----------------" << std::endl;
    for(unsigned int i=0; i<KnotVectorX.size(); i++){
        std::cout << KnotVectorX[i] << ", ";
    }
    std::cout << std::endl;

    std::cout<< "---------------IntervalsY-----------------" << std::endl;
    
    for(unsigned int i=0; i<IntervalsY.size(); i++){
        std::cout << IntervalsY[i] << ", ";
    }
    std::cout << std::endl;

    std::cout<< "---------------IntervalsX-----------------" << std::endl;
    
    for(unsigned int i=0; i<IntervalsX.size(); i++){
        std::cout << IntervalsX[i] << ", ";
    }
    std::cout << std::endl;
}


//int nBasisX=NewCtrlPtsW.size();
//int nBasisY=NewCtrlPtsW[0].GetNumCols();
int nX= KnotVectorX.size()-pX-2;
int nY= KnotVectorY.size()-pY-2;
int nLoc=(pX+1)*(pY+1);
std::vector<Eigen::Triplet<double>> tripletList;
tripletList.reserve(hX * hY * nLoc * nLoc);

//for paralelitation
int nthreads = omp_get_max_threads();
std::vector<std::vector<Eigen::Triplet<double>>> threadTriplets(nthreads);
std::vector<std::vector<double>> threadLinearForms(nthreads, std::vector<double>(size, 0.0));

// Heuristica de reserva por hilo para evitar muchas realocaciones.
// Ajusta segun memoria / experiencia.
for (int t = 0; t < nthreads; ++t) {
    threadTriplets[t].reserve((hX * hY * nLoc * nLoc) / nthreads / 4 + 1000);
}

Eigen::MatrixXd GaussWeights = weightsX * weightsY.transpose();
/*
cout << "weightsY: " <<weightsY << std::endl;
cout << "weightsX: " << weightsX << std::endl;
cout << "GaussWeights: " << std::endl << GaussWeights << std::endl;
*/
auto start = std::chrono::high_resolution_clock::now();
//Precomputation of BasisFuns and its derivates.
std::vector<int> SpanIndexY(IntervalsY.size()-1);
std::vector<int> SpanIndexX(IntervalsX.size()-1);

std::vector<std::vector<iMatrix<double>>> Basis_and_DersY=PreBasis_and_Ders(pY, 1, KnotVectorY, IntervalsY, nodesY,SpanIndexY);
std::vector<std::vector<iMatrix<double>>> Basis_and_DersX=PreBasis_and_Ders(pX, 1, KnotVectorX, IntervalsX, nodesX,SpanIndexX);

if(ShowValues){
    std::cout<< "---------------Basis_and_DersY-----------------" << std::endl;
    for(unsigned int i=0; i< Basis_and_DersY.size(); i++){
        std::cout<< "-   -   -   -   -   -   -   -   -   -   -   -" << std::endl;
        for(unsigned int j=0; j< nEvalsY; j++){
            Basis_and_DersY[i][j].PrintMatrix();
            std::cout<< "-   -   -   -   -   -   -   -   -   -   -   -" << std::endl;
        }
        std::cout<< "-+--+--+--+--+--+--+--+--+--+--+-" << std::endl;
    }
    std::cout<< "---------------Basis_and_DersX-----------------" << std::endl;
    for(unsigned int i=0; i< Basis_and_DersX.size(); i++){
        std::cout<< "-   -   -   -   -   -   -   -   -   -   -   -" << std::endl;
        for(unsigned int j=0; j< nEvalsX; j++){
            Basis_and_DersX[i][j].PrintMatrix();
            std::cout<< "-   -   -   -   -   -   -   -   -   -   -   -" << std::endl;
        }
        std::cout<< "-+--+--+--+--+--+--+--+--+--+--+-" << std::endl;
    }

}
//std::unordered_map<int, std::vector<std::pair<double,double>>> spanMapY=BuildSpanMap(IntervalsY, SpanIndexY);
//std::unordered_map<int, std::vector<std::pair<double,double>>> spanMapX=BuildSpanMap(IntervalsX, SpanIndexX);

// Double index matrix for easier and faster assembly
std::vector<double> detg_vals;
#pragma omp parallel for collapse(2) schedule(static)
// Assembly of the matrix
for (unsigned int i = 0; i < hX; i++) {
    for (unsigned int j = 0; j < hY; j++) {

        int tid = omp_get_thread_num();
        // DEBUG: vector local por hilo
        std::vector<double> local_detg_vals;


        // Local buffers referidos por tid
        auto &localTriplets = threadTriplets[tid];
        auto &localLinear   = threadLinearForms[tid];

        iMatrix<double> K_e(subMatSize, subMatSize);
        std::vector<double> F_e((pX+1)*(pY+1));
        if (ShowIterations) {
            std::cout<<"--------Active Intervals---------"<< std::endl;
            std::cout<< "IntervalX: (" << IntervalsX[i] << " " << IntervalsX[i+1]<<")"<< std::endl;
            std::cout<< "IntervalY: (" << IntervalsY[j] << " " << IntervalsY[j+1]<<")"<< std::endl;
            std::cout<< "SpanIndexX: (" << SpanIndexX[i] <<")"<< std::endl;
            std::cout<< "SpanIndexY: (" << SpanIndexY[j] <<")"<< std::endl;
        }
        double auxeval1X= (IntervalsX[i+1] - IntervalsX[i])/2;
        double auxeval2X= (IntervalsX[i+1] + IntervalsX[i])/2;
        double auxeval1Y= (IntervalsY[j+1] - IntervalsY[j])/2;
        double auxeval2Y= (IntervalsY[j+1] + IntervalsY[j])/2;

        
        int spanU = SpanIndexX[i];
        int spanV = SpanIndexY[j];
        iMatrix<double> activeWeights=NewWeights.GetSubMat(spanU-pX,spanU, spanV-pY, spanV).Transpose();
        auto &DersBasisX_vec = Basis_and_DersX[i]; // vector<iMatrix<double>> size nEvalsX
        auto &DersBasisY_vec = Basis_and_DersY[j]; // vector<iMatrix<double>> size nEvalsY
        //std::vector<double> WeightsY=NewWeights.GetSubMat(0,NewWeights.GetNumRows()-1, spanV-pY, spanV-pY).GetCol(0);
        //std::vector<double> WeightsX=NewWeights.GetSubMat(spanU-pX, spanU-pX,0,NewWeights.GetNumCols()-1).GetRow(0);
        //std::vector<iMatrix<double>> RationalizedDersBasisY=RationalizeDersBasisFunsVec(spanV, nEvalsY, nodesY, IntervalsY[j],IntervalsY[j+1], pY, 1, nY,
        // KnotVectorY, Basis_and_DersY[j], WeightsY);
        //std::vector<iMatrix<double>> RationalizedDersBasisX=RationalizeDersBasisFunsVec(spanU, nEvalsX,nodesX, IntervalsX[i], IntervalsX[i+1],pX, 1, nX ,
        // KnotVectorX, Basis_and_DersX[i], WeightsX);
        if (ShowIterations2) {
            std::cout<<"--------RationalizedDersBasisY---------"<< std::endl;
            for(unsigned int iter=0; iter<nEvalsY; iter++){
                //RationalizedDersBasisY[iter].PrintMatrix();
                std::cout<<"------------------------------"<< std::endl;
            }
            std::cout<<"--------RationalizedDersBasisX---------"<< std::endl;
            for(unsigned int iter=0; iter<nEvalsX; iter++){
                //RationalizedDersBasisX[iter].PrintMatrix();
                std::cout<<"------------------------------"<< std::endl;
            }
        }

        //Emulation of computations of the matrix
        //iMatrix<iMatrix<double>> SurfBasisRat, SurfBasisRatX, SurfBasisRatY;
        for(int a = 0; a < nEvalsX; ++a){

            double EvalPointX=auxeval1X*nodesX(a)+auxeval2X;
            //std::cout << "EvalPointX: " << EvalPointX << ". ";
            for(int b = 0; b < nEvalsY; ++b){
                double EvalPointY=auxeval1Y*nodesY(b)+auxeval2Y;
                //std::cout << "EvalPointY: " << EvalPointY << ". ";
                iMatrix<std::vector<double>> Allders = SurfaceDerivsAlg1(KnotVectorY.size()-2-pY,pY,KnotVectorY,KnotVectorX.size()-2-pX,pX,KnotVectorX,NewCtrlPtsW,
                                                                         EvalPointY, EvalPointX, 1);
                
                
                iMatrix<std::vector<double>> Aders(2, 2);
                iMatrix<double> Wders(2, 2);
                for(unsigned int iteri=0; iteri<=1; iteri++){
                    for(unsigned int iterj=0; iterj<=1; iterj++){
                        Aders(iteri,iterj)= std::vector<double>(Allders(iteri,iterj).begin(), Allders(iteri,iterj).begin()+3);
                        Wders(iteri,iterj)=Allders(iteri,iterj)[3];
                    }
                }
                iMatrix<iMatrix<double>> RationalizedBasis2d=RationalizeDersBasisFuns2D(DersBasisX_vec[a], DersBasisY_vec[b],activeWeights,pX,pY);
                std::vector<double> SurfBasisRat=(RationalizedBasis2d(0,0)).Vectorize();
                std::vector<double> SurfBasisRatX=(RationalizedBasis2d(1,0)).Vectorize();
                std::vector<double> SurfBasisRatY=(RationalizedBasis2d(0,1)).Vectorize();
                iMatrix<std::vector<double>> SurfDers= RatSurfaceDerivs(Aders, Wders, 1);
                
                iMatrix<double> PartialXX=SurfBasisRatX*SurfBasisRatX;
                iMatrix<double> PartialXY=SurfBasisRatX*SurfBasisRatY;
                iMatrix<double> PartialYY=SurfBasisRatY*SurfBasisRatY;
                double g11=dot_product(SurfDers(1,0),SurfDers(1,0));
                double g12=dot_product(SurfDers(1,0),SurfDers(0,1));
                double g22=dot_product(SurfDers(0,1),SurfDers(0,1));
                double detg=g11*g22-g12*g12;
                double J=std::sqrt(detg);
                double Jinv=1/J;
                double dOmega = auxeval1X * auxeval1Y * GaussWeights(a,b);
                double dA = dOmega*J;
                double auxconstant = dOmega * Jinv;
                double g11c=  g22*auxconstant;
                double g12c= -g12*auxconstant;
                double g22c=  g11*auxconstant;
                //std::cout << "<<<<<<<"<< std::endl;
                local_detg_vals.push_back(detg);
                //iMatrix<double> Kq = g22c*PartialXX + g12c*(PartialXY + PartialXY.Transpose()) + g11c*PartialYY;
                iMatrix<double> Kq(subMatSize, subMatSize);
                Kq = PartialXX;
                Kq *= g22c;

                iMatrix<double> Temp = PartialXY + PartialXY.Transpose(); 
                Temp *= g12c; 
                Kq += Temp;

                Temp = PartialYY; 
                Temp *= g11c; 
                Kq += Temp; 

                K_e +=  Kq;
                for (int A = 0; A < SurfBasisRat.size(); ++A) {
                    F_e[A] += SurfBasisRat[A] * f(EvalPointX, EvalPointY) * dA;
                }
            }
        }

        // Emulation of local load vector (all ones)
        
        // Rellenar triplets en matriz global
        for(int a = 0; a <= pX; ++a){
            
            for(int b = 0; b <= pY; ++b){
                int i_local = b*(pX+1) + a; 
                int I_global = (spanU - pX + a) + (spanV - pY + b)*nElementsX;
                localLinear[I_global] += F_e[i_local];
                //std::cout<< "---------------------------"<< std::endl;
                //std::cout << "I_global= " << I_global << std::endl;

                //LinearForm[I_global] += F_e[i_local];

                for(int c = 0; c <= pX; ++c){
                    for(int d = 0; d <= pY; ++d){
                        
                        int j_local = d*(pX+1) + c; 
                        int J_global = (spanU - pX + c) + (spanV - pY + d)*nElementsX;
                        
                        //std::cout << "J_global= " << J_global << std::endl;
                        double val = K_e(i_local, j_local);
                        if(val != 0.0){
                            //localTriplets.emplace_back(I_global, J_global, val);
                            localTriplets.emplace_back(I_global, J_global, val);
                        }
                    }
                }
            }
        }

        #pragma omp critical
        {
            detg_vals.insert(detg_vals.end(),
                                local_detg_vals.begin(),
                                local_detg_vals.end());
        }
    }
}

for (int t = 0; t < nthreads; ++t) {
    // mover triplets al tripletList (efficient push_back many)
    tripletList.insert(tripletList.end(), threadTriplets[t].begin(), threadTriplets[t].end());

    // sumar fuerza local en LinearForm
    for (int idx = 0; idx < size; ++idx) {
        LinearForm[idx] += threadLinearForms[t][idx];
    }
}
global.setFromTriplets(tripletList.begin(), tripletList.end());

auto end = std::chrono::high_resolution_clock::now();
auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
std::cout << "Tiempo tomado: " << duration.count() << " nanosegundos\n";

if(PrintMatrixes){
    std::cout<< "--------------GlobalMat----------------" << std::endl;
    for (int k = 0; k < global.outerSize(); ++k){
        for (Eigen::SparseMatrix<double>::InnerIterator it(global, k); it; ++it){
            std::cout << "(" << it.row() << "," << it.col() << "): " << it.value() << "\n";
        }
    }
    //std::cout<<Eigen::MatrixXd(global) <<std::endl;
}

if(PrintMatrixes){
    std::cout<< "--------------GlobalLinearForm----------------" << std::endl;
    for (int k = 0; k < LinearForm.size(); ++k){
        std::cout   << LinearForm[k] << ", ";
    }
    std::cout << std::endl;
}

// -------------------- Imposicion de Dirichlet homogeneo (u = 0) --------------------
std::vector<char> isDirichlet(size, false);
std::vector<int> dirichletNodes;
dirichletNodes.reserve(2 * (nElementsX + nElementsY));

for (int j = 0; j < nElementsY; ++j) {
    for (int i = 0; i < nElementsX; ++i) {
        int idx = i + j * nElementsX;
        if (i == 0 || i == nElementsX - 1 || j == 0 || j == nElementsY - 1) {
            isDirichlet[idx] = 1;
            dirichletNodes.push_back(idx);
        }
    }
}
if(ShowConditions){
    std::cout << "Numero de nodos Dirichlet detectados: " << dirichletNodes.size() << std::endl;
}

for (int col = 0; col < global.outerSize(); ++col) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(global, col); it; ++it) {
        int row = it.row();
        int c   = it.col();
        if (isDirichlet[row] && row != c) {
            it.valueRef() = 0.0;
        }
    }
}

for (int k = 0; k < (int)dirichletNodes.size(); ++k) {
    int colIdx = dirichletNodes[k];
    for (Eigen::SparseMatrix<double>::InnerIterator it(global, colIdx); it; ++it) {
        int row = it.row();
        if (row != colIdx) {
            it.valueRef() = 0.0;
        }
    }
}

for (int k = 0; k < (int)dirichletNodes.size(); ++k) {
    int idx = dirichletNodes[k];
    global.coeffRef(idx, idx) = 1.0;
    LinearForm[idx] = 0.0;
}

if(PrintMatrixes2){
    
    std::cout<< "--------------GlobalMat----------------" << std::endl;
    for (int k = 0; k < global.outerSize(); ++k){
        for (Eigen::SparseMatrix<double>::InnerIterator it(global, k); it; ++it){
            std::cout << "(" << it.row() << "," << it.col() << "): " << it.value() << "\n";
        }
    }
    std::cout<< "--------------GlobalLinearForm----------------" << std::endl;
    for (int k = 0; k < LinearForm.size(); ++k){
        std::cout   << LinearForm[k] << ", ";
    }
    std::cout << std::endl;

}

global.makeCompressed();

if(ShowSolutions){
    std::cout << "Resolviendo sistema (K u = F) con Conjugate Gradient..." << std::endl;
}

int maxIter = std::max(100, std::min(5 * size, 500000));
auto tsolve_start = std::chrono::high_resolution_clock::now();
Eigen::VectorXd solution;
// IncompleteCholesky como precondicionador (SPD)
try {
    using Precond = Eigen::IncompleteCholesky<double>;
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Precond> solver;
    solver.compute(global);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Precond compute fallo (IC).");
    }
    solver.setTolerance(1e-8);
    solver.setMaxIterations(maxIter);
    solution = solver.solve(LinearForm);
    auto tsolve_end = std::chrono::high_resolution_clock::now();

    double elapsed_s = std::chrono::duration_cast<std::chrono::duration<double>>(tsolve_end - tsolve_start).count();
    std::cout << "Tiempo solucion: " << elapsed_s << " s\n";
    if(ShowSolutions){
        std::cout << "Solver info (IC precond): " << solver.info() << "\n";
        std::cout << "Iterations: " << solver.iterations() << ", Estimated error: " << solver.error() << "\n";
        std::cout << "Solucion: " << solution.transpose() << std::endl;
    }
    
    

} catch (const std::exception &e) {
    // Fallback: usar precondicionador diagonal (Jacobi) si IC falla
    if(ShowSolutions){
    std::cerr << "Fallo precondicionador IncompleteCholesky o excepcion: " << e.what() << "\n";
    std::cerr << "Usando ConjugateGradient + DiagonalPreconditioner (fallback)." << std::endl;
    }
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::DiagonalPreconditioner<double>> solver2;
    solver2.compute(global);
    solver2.setTolerance(1e-8);
    solver2.setMaxIterations(2*maxIter);
    solution = solver2.solve(LinearForm);

    auto tsolve_end2 = std::chrono::high_resolution_clock::now();
    double elapsed_s2 = std::chrono::duration_cast<std::chrono::duration<double>>(tsolve_end2 - tsolve_start).count();
    std::cout << "Tiempo solucion: " << elapsed_s2 << " s\n";
    if(ShowSolutions){
        std::cout << "Solver info (Diag precond): " << solver2.info() << "\n";
        std::cout << "Iterations: " << solver2.iterations() << ", Estimated error: " << solver2.error() << "\n";
        std::cout << "Solucion: " << solution.transpose() << std::endl;
    }
}

// -------------------- Evaluacion en los elementos de la Solucion --------------------

    
        
iMatrix<double> MatSol(nElementsX, nElementsY, solution);
//iMatrix<double> MatSol=MatSol1.Transpose();
iMatrix<double> SolEvals(nEvalsX*hX, nEvalsY*hY);
iMatrix<double> SolEvalsX(nEvalsX*hX, nEvalsY*hY);
iMatrix<double> SolEvalsY(nEvalsX*hX, nEvalsY*hY);
iMatrix<double> SolAnalytic(nEvalsX*hX, nEvalsY*hY);
iMatrix<double> SolAnalyticX(nEvalsX*hX, nEvalsY*hY);
iMatrix<double> SolAnalyticY(nEvalsX*hX, nEvalsY*hY);
iMatrix<std::vector<double>> Surf(nEvalsX*hX, nEvalsY*hY);
for (unsigned int i = 0; i < hX; i++) {
    for (unsigned int j = 0; j < hY; j++) {

        double auxeval1X= (IntervalsX[i+1] - IntervalsX[i])/2;
        double auxeval2X= (IntervalsX[i+1] + IntervalsX[i])/2;
        double auxeval1Y= (IntervalsY[j+1] - IntervalsY[j])/2;
        double auxeval2Y= (IntervalsY[j+1] + IntervalsY[j])/2;

        
        int spanU = SpanIndexX[i];
        int spanV = SpanIndexY[j];
        iMatrix<double> activeWeights=NewWeights.GetSubMat(spanU-pX,spanU, spanV-pY, spanV).Transpose();
        auto &DersBasisX_vec = Basis_and_DersX[i]; // vector<iMatrix<double>> size nEvalsX
        auto &DersBasisY_vec = Basis_and_DersY[j]; // vector<iMatrix<double>> size nEvalsY
        for(int a = 0; a < nEvalsX; ++a){

            double EvalPointX=auxeval1X*nodesX(a)+auxeval2X;
            //std::cout << "EvalPointX: " << EvalPointX << ". ";
            for(int b = 0; b < nEvalsY; ++b){
                double EvalPointY=auxeval1Y*nodesY(b)+auxeval2Y;
                iMatrix<iMatrix<double>> RationalizedBasis2d=RationalizeDersBasisFuns2D(DersBasisX_vec[a], DersBasisY_vec[b],activeWeights,pX,pY);
                
                iMatrix<double> Rbasis = RationalizedBasis2d(0, 0);
                iMatrix<double> RbasisX = RationalizedBasis2d(1, 0);
                iMatrix<double> RbasisY = RationalizedBasis2d(0, 1);
                iMatrix<double> LocalSol = MatSol.GetSubMat(spanU - pX, spanU, spanV - pY, spanV);
                double u_h = 0.0;
                double u_hX= 0.0;
                double u_hY= 0.0;
                for (int ii = 0; ii < Rbasis.GetNumRows(); ++ii) {
                    for (int jj = 0; jj < Rbasis.GetNumCols(); ++jj) {
                        u_h += Rbasis(ii, jj) * LocalSol(jj, ii);
                        u_hX+= RbasisX(ii,jj) * LocalSol(jj, ii);
                        u_hY+= RbasisY(ii,jj) * LocalSol(jj, ii);
                    }
                }
                SolEvals(i * nEvalsX + a, j * nEvalsY + b) = u_h;
                SolEvalsX(i * nEvalsX + a, j * nEvalsY + b) = u_hX;
                SolEvalsY(i * nEvalsX + a, j * nEvalsY + b) = u_hY;
                std::vector<double> SurfPoint=SurfacePointRational(nY, pY, KnotVectorY, nX, pX, KnotVectorX, NewCtrlPtsW, EvalPointY, EvalPointX);
                Surf(i*nEvalsX + a, j*nEvalsY + b)=SurfPoint;
                SolAnalytic(i*nEvalsX+a, j*nEvalsY +b) = analyticSol(SurfPoint[0], SurfPoint[1]);
                SolAnalyticX(i*nEvalsX+a, j*nEvalsY +b) = analyticSol_dx(SurfPoint[0], SurfPoint[1]);
                SolAnalyticY(i*nEvalsX+a, j*nEvalsY +b) = analyticSol_dy(SurfPoint[0], SurfPoint[1]);
                /*
                std::cout<<"----------------------------------"<<std::endl;
                std::cout<< "SurfPoint[0]=" << SurfPoint[0]<< ", EvalPointX="<< EvalPointX <<std::endl;
                std::cout<< "SurfPoint[1]=" << SurfPoint[1]<< ", EvalPointY="<< EvalPointY <<std::endl;
                std::cout << "Analytic= " << SolAnalytic(i*nEvalsX+a, j*nEvalsY +b) << ", Numeric= " << SolEvals(i*nEvalsX + a, j*nEvalsY + b) << std::endl;
                */
            }
        }
    }
}
//MatSol.PrintMatrix();

//ExportToTxt(SolEvals, Surf, PATH + "surfaceCircle3d_hX="+  std::to_string(static_cast<int>(hX)) + "_pX=" + to_string(pX) +"_hY="+  std::to_string(static_cast<int>(hY)) +
//                               "_p=" + to_string(pY) + " data.txt");

//------------------L2 norm calculations----------------------

double L2_error_sq = 0.0;
double H1_semi_error_sq = 0.0;

for (unsigned int i = 0; i < hX; i++) {
    double auxeval1X= (IntervalsX[i+1] - IntervalsX[i])/2;
    double auxeval2X= (IntervalsX[i+1] + IntervalsX[i])/2;
    for (unsigned int j = 0; j < hY; j++) {
        double auxeval1Y= (IntervalsY[j+1] - IntervalsY[j])/2;
        double auxeval2Y= (IntervalsY[j+1] + IntervalsY[j])/2;
        double h_elementos = auxeval1X*auxeval1Y;
        double cumSumElements=0.0;
        double cumSumH1 = 0.0;
        for (int a = 0; a < nEvalsX; ++a) {
            double EvalPointX=auxeval1X*nodesX(a)+auxeval2X;
            for (int b = 0; b < nEvalsY; ++b) {
                // --- Diferencial parametrico ---

                
            
            
                double EvalPointY=auxeval1Y*nodesY(b)+auxeval2Y;
                iMatrix<std::vector<double>> Allders = SurfaceDerivsAlg1(KnotVectorY.size()-2-pY,pY,KnotVectorY,KnotVectorX.size()-2-pX,pX,
                                                                        KnotVectorX,NewCtrlPtsW, EvalPointY, EvalPointX, 1);
                
                
                iMatrix<std::vector<double>> Aders(2, 2);
                iMatrix<double> Wders(2, 2);
                for(unsigned int iteri=0; iteri<=1; iteri++){
                    for(unsigned int iterj=0; iterj<=1; iterj++){
                        Aders(iteri,iterj)= std::vector<double>(Allders(iteri,iterj).begin(), Allders(iteri,iterj).begin()+3);
                        Wders(iteri,iterj)=Allders(iteri,iterj)[3];
                    }
                }
                iMatrix<std::vector<double>> SurfDers= RatSurfaceDerivs(Aders, Wders, 1);
                
                // Punto fisico (x, y, z) de la superficie S
                std::vector<double> SurfPoint = SurfDers(0, 0);
                double x_phys = SurfPoint[0];
                double y_phys = SurfPoint[1];

                // Gradiente fisico de la solucion analitica (grad x u)
                double du_dx_analytic = -x_phys / 2.0;
                double du_dy_analytic = -y_phys / 2.0;

                // Gradiente parametrico de la solucion analitica (grad xi u) usando la Regla de la Cadena
                // dS/dxi = SurfDers(0, 1) = (dx/dxi, dy/dxi, dz/dxi)
                std::vector<double> S_xi = SurfDers(0, 1);
                double du_dxi_analytic = du_dx_analytic * S_xi[0] + du_dy_analytic * S_xi[1];
                // (El termino dz/dxi * du/dz es 0 ya que u no depende de z)

                // dS/deta = SurfDers(1, 0) = (dx/deta, dy/deta, dz/deta)
                std::vector<double> S_eta = SurfDers(1, 0);
                double du_deta_analytic = du_dx_analytic * S_eta[0] + du_dy_analytic * S_eta[1];


                double g11=dot_product(SurfDers(1,0),SurfDers(1,0));
                double g12=dot_product(SurfDers(1,0),SurfDers(0,1));
                double g22=dot_product(SurfDers(0,1),SurfDers(0,1));
                
                double detg=g11*g22-g12*g12;
                double J=std::sqrt(detg);
                double Jinv=1/J;

                double g11c_J = g22 * Jinv; // g^11 * J
                double g12c_J = -g12 * Jinv; // g^12 * J
                double g22c_J = g11 * Jinv; // g^22 * J
                //double dXiEta = auxeval1X * auxeval1Y * GaussWeights(a,b);
                double diff = SolEvals(i*nEvalsX + a, j*nEvalsY + b) - SolAnalytic(i*nEvalsX + a, j*nEvalsY + b);
                double diff_dx = SolEvalsX(i*nEvalsX + a, j*nEvalsY + b) - du_dxi_analytic; // d/dxi
                double diff_dy = SolEvalsY(i*nEvalsX + a, j*nEvalsY + b) - du_deta_analytic; // d/deta
                

                // --- Aportacion al error L2 ---
                cumSumElements += diff * diff * GaussWeights(a,b)*J;
                cumSumH1 += (diff_dx * diff_dx *g22c_J +2.0*g12c_J*diff_dx*diff_dy + diff_dy * diff_dy*g11c_J) * GaussWeights(a,b);
            }
        }
        L2_error_sq+=cumSumElements*h_elementos;
        H1_semi_error_sq += cumSumH1 * h_elementos;
    }
}

double L2_error = std::sqrt(L2_error_sq);
double H1_error  = sqrt(L2_error_sq + H1_semi_error_sq);
std::cout << "Norma L2 = " << L2_error << std::endl;
std::cout << "Norma H1 = " << H1_error << std::endl;