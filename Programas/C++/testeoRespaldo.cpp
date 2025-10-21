/******************************************************************************* 
 *AUTHOR: Carlos Palomera Oliva
 *DATE: MARCH 2025
 *******************************************************************************/

// C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\Programas\C++
// compilar con avx2 -mavx2
// g++ testeo.cpp iMatrix.cpp funciones.cpp SparseTensor.cpp -fopenmp -O3 -march=native -o testeo
//#include"SparseTensor.cpp" 
#include <iostream>
#include <omp.h>
#include <vector>
#include"iMatrixSecurity.cpp"
#include"funciones.hpp"
#include <fstream>
#include <chrono>
#include <cstdint>
using namespace std;
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;


inline uint64_t rdtsc() {
    unsigned int lo, hi;
    asm volatile (
        "rdtsc" 
        : "=a" (lo), "=d" (hi)
    );
    return ((uint64_t)hi << 32) | lo;
}

int main(int argc, char *argv[]) {
    


    if(stoi(argv[1])==1){//test superficies de Bezier
    /*
    iMatrix matrix = ReadDataFile_Matrix("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/Programas/C++/EjemploFicheroPtos2d.txt");
    matrix.PrintMatrix();

    matrix=BezierCurve(matrix, 10);
    matrix.PrintMatrix();
    /************************************************************ */
    /*
    iMatrix matrix2 = ReadDataFile_Matrix("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/Programas/C++/EjemploFicheroPtos3d.txt");
    matrix2.PrintMatrix();

    matrix2=BezierCurve(matrix2, 10);
    matrix2.PrintMatrix();
    */
   std::vector<iMatrix<double>> CtrlPts = ReadDataFile_CtrlPts("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/EjemploFicheroTensor.txt");
   /*
   for(unsigned int i=0; i<Tensor.size();i++){
       Tensor[i].PrintMatrix();
       std::cout<<"-------------------" << std::endl;
   }
   */
   /*
   Tensor=BezierSurface(Tensor, 10, 10);
   for(unsigned int i=0; i<10;i++){
       Tensor[i].PrintMatrix();
   }
   */  
   //Tensor = ReadDataFile_Tensor("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/EjemploFicheroTensor.txt");
   std::vector<iMatrix<double>> Tensor = BezierSurface(CtrlPts, 1000);

   /*
   for(unsigned int i=0; i<100;i++){
       Tensor[i].PrintMatrix();
   }
   */
   Write_BezierSurface(Tensor, "C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/BezierSurf/TestFicheroDatosBezierSurface");
   Write_CtrlPts(CtrlPts,"C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/BezierSurf/TestFicheroDatosPtosCtrl");
   std::cout<<"fin 1"<<endl;
   /*
   iMatrix<double> matriz2(3,4);
   matriz2.PrintMatrix();
   */


   /*
   //Test sobre la clase matriz
   iMatrix<double> matriz;
   matriz.Resize(2,2);
   matriz.SetToIdentity();
   matriz.PrintMatrix();

   bool x = matriz(1,1)= 37);
   x = matiz(1,0)= 14);
   matriz(,1)= 7);
   matriz.PintMatrix();

   cout<<"llega"<<endl;
   cout << "columna0 "<< matriz.GetCol(0)[0] << " " << matriz.GetCol(0)[1] <<endl;
   cout << "columna1 "<< matriz.GetCol(1)[0] << " " << matriz.GetCol(1)[1] <<endl;    

   cout << "fila0 "<< matriz.GetRow(0)[0] << " " << matriz.GetRow(0)[1] <<endl;
   cout << "fila1 "<< matriz.GetRow(1)[0] << " " << matriz.GetRow(1)[1] <<endl;    
   */
    }

    else if(stoi(argv[1])==3){ //test find span y BasisFuns

        int p=3;
        std::vector<double> KnotVector={0,0,0,0.2,0.4,0.6,0.8,1,1,1};
        /*
        std::vector<double> t_values ={0.1,0.54,1.35,1.79,2,3,4,4.5,5,8};
        for(unsigned int iter=0; iter<t_values.size();iter++){
            auto it = std::upper_bound(Nurbs.begin(), Nurbs.end(), t_values[iter]); //componente no nula
            int itaux =static_cast<int>(std::distance(Nurbs.begin(), it))-1;
            int i = FindSpan(Nurbs,t_values[iter], p, Nurbs.size()-2-p);
            std::cout << t_values[iter]<<" t: " << i<< " i"<<std::endl;
            std::cout << t_values[iter]<<" t: " << itaux<< " iter"<<std::endl;
            
        }
        */
        double N_Steps = 100;

        iMatrix<double> sol(KnotVector.size()-p,N_Steps+1);
        iMatrix<double> sol2(KnotVector.size()-p,N_Steps+1);
        //std::vector<double> Nurb(NurbsVector.size());
        double Step=1/(N_Steps);
        for(unsigned int i = 0; i<N_Steps+1; i++){
            double t=Step*i*8;
            //std::cout << "Iteración iniciada " << i<< std::endl;
            int span = FindSpan(KnotVector, t, p, KnotVector.size()-2-p);
            std::vector<double> eval= BasisFuns(span, t , p, KnotVector);
            
            for(int j = span-p; j<=span;j++){
                if(j>=0) {
                    sol(j,i)=eval[j-span+p];                }
                else{
                    std::cout<<"Función base con indices " << j<<" "<< p << " no ha sido dibujada en el punto t="<< t << std::endl;
                    std::cout<<"Valor de la función previa: "<< eval[j-span+p]<< std::endl;
                }
                //std::cout<< "j: " << j<< " span " << span << " t:" << t << " i: " << i<<std::endl;
                
            }
            std::cout<< " span " << span << " t:" << t <<std::endl;
            iMatrix<double> AllBasis=AllBasisFuns(span, t , p, KnotVector);
            AllBasis.PrintMatrix();
            sol2.SetCol(i,AllBasisFuns(span, Step*i , p, KnotVector).GetRow(p));
            //std::cout << "Iteración completada " << i<< std::endl;
        }
        Write_BasisFunctions(sol, "C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/BasisFuncts");
        Write_BasisFunctions(sol2, "C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/AllBasisFuncts");




        std::cout<<"fin 3"<<endl;


    }


    else if(stoi(argv[1])==4){//test DerBasisFuns
        int p=3;
        std::vector<double> Nurbs={0,0,0,0,2,4,4,6,6,6,8,8,8,8};


        double N_Steps = 100;
        std::vector<double> eval(Nurbs.size()-p);
        iMatrix<double> ders(p+1,p+1);
        iMatrix<double> sol(Nurbs.size()-p,N_Steps+1);
        //std::vector<double> Nurb(NurbsVector.size());
        double Step=1/(N_Steps);
        auto start = std::chrono::high_resolution_clock::now();
        for(unsigned int i = 0; i<N_Steps+1; i++){
            
            double t=Step*i*8;
            //std::cout << "Iteración iniciada " << i<< std::endl;
            int span = FindSpan(Nurbs, t, p, Nurbs.size()-2-p);
            ders = DerBasisFuns(span,t,p,1,Nurbs);
            ders.PrintMatrix();
            std::cout<<"---------------------"<<std::endl;
            for(int j = span-p; j<=span;j++){
                if(j>=0) {
                    sol(j,i)=ders(1,j-span+p);
                }
                else{
                    std::cout<<"Función base con indices " << j<<" "<< p << " no ha sido dibujada en el punto t="<< t << std::endl;
                    std::cout<<"Valor de la función previa: "<< eval[j-span+p]<< std::endl;
                }
                //std::cout<< "j: " << j<< " span " << span << " t:" << t << " i: " << i<<std::endl;
                
            }
        }
        
        Write_BasisFunctions(sol, "C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/BasisFuncts_Derivates");
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    
        std::cout << "Tiempo tomado: " << duration.count() << " nanosegundos\n";
    
        std::cout<<"fin 4"<<endl;
    
    }

    else if(stoi(argv[1])==5){//test DerBasisFuns
        int p=3;
        std::vector<double> Nurbs={0,0,0,0,2,4,4,6,6,6,8,8,8,8};


        double N_Steps = 10000;
        std::vector<double> eval(Nurbs.size()-p);
        iMatrix<double> ders(p+1,p+1);
        iMatrix<double> sol(Nurbs.size()-p,N_Steps+1);
        //std::vector<double> Nurb(NurbsVector.size());
        double Step=1/(N_Steps);
        uint64_t start = rdtsc();
        for(unsigned int i = 0; i<N_Steps+1; i++){
            
            double t=Step*i*8;
            //std::cout << "Iteración iniciada " << i<< std::endl;
            int span = FindSpan(Nurbs, t, p, Nurbs.size()-2-p);
            ders = DerBasisFuns(span,t,p,1,Nurbs);
            //ders.PrintMatrix();
            //std::cout<<"---------------------"<<std::endl;
            for(int j = span-p; j<=span;j++){
                if(j>=0) {
                    sol(j,i)=ders(1,j-span+p);
                }
                else{
                    std::cout<<"Función base con indices " << j<<" "<< p << " no ha sido dibujada en el punto t="<< t << std::endl;
                    std::cout<<"Valor de la función previa: "<< eval[j-span+p]<< std::endl;
                }
                //std::cout<< "j: " << j<< " span " << span << " t:" << t << " i: " << i<<std::endl;
                
            }
        }
        
        Write_BasisFunctions(sol, "C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/BasisFuncts_Derivates");
        uint64_t end = rdtsc();
        uint64_t cycles = end - start;

        std::cout << "Ciclos de reloj tomados: " << cycles << "\n";
    
        std::cout<<"fin 5"<<endl;
    
    }

    else if(stoi(argv[1])==6){//test CurvePoint
        
        iMatrix matrix = ReadDataFile_Matrix("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/Ej1.txt");
        //cout<<"abierto"<<endl;
        matrix.PrintMatrix();
        std::vector<double> KnotVector={0,0,0,0.2,0.4,0.6,0.8,1,1,1};
        int p=2;
        int dim=2;
        double N_Steps = 100;

        iMatrix<double> sol(dim,N_Steps+1);
        std::vector<double> eval(dim);
        double Step=1/(N_Steps);
        for(unsigned int i = 0; i<N_Steps+1; i++){
            
            eval=CurvePoint(KnotVector.size()-2-p,p,KnotVector,matrix,i*Step);
            sol.SetCol(i,eval);
        }
        Write_Curve(sol,"C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/Curvas/CurvaEj1.txt");
        std::cout<<"fin 6"<<endl;
    }
    
    
    else if(stoi(argv[1])==7){//test CurveDerivsAlg1
        iMatrix matrix = ReadDataFile_Matrix("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/Ej1.txt");
        //cout<<"abierto"<<endl;
        matrix.PrintMatrix();
        std::vector<double> Nurbs={0,0,0,0,0.25,0.5,0.75,1,1,1,1};
        int p=3;
        int dim=2;
        double N_Steps = 100000;

        std::vector<iMatrix<double>> sol(N_Steps+1);
        iMatrix<double> eval(matrix.GetNumRows(),std::min(dim,p)+1);
        double Step=1/(N_Steps);
        auto start = std::chrono::high_resolution_clock::now();
        for(unsigned int i = 0; i<N_Steps+1; i++){
            eval=CurveDerivsAlg1(Nurbs.size()-2-p,p,Nurbs,matrix, i*Step,dim-1);
            //eval=CurveDerivsAlg1(Nurbs.size()-2-p,p,Nurbs,matrix,i*Step);
            sol[i]=eval;
            //sol[i].PrintMatrix();
            //std::cout<< "----------------" <<std::endl;
            //eval.PrintMatrix();
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        iMatrix<double> Derivadaiesima(2, N_Steps+1);
        for(unsigned int iter=0; iter<N_Steps+1;iter++){
            Derivadaiesima.SetCol(iter,sol[iter].GetRow(1));
        }
        std::cout << "Tiempo tomado: " << duration.count() << " nanosegundos\n";
        Write_Curve(Derivadaiesima ,"C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/Curvas/CurvaEj1DersAlg1.txt");
        std::cout<<"fin 7"<<endl;
    }

    else if(stoi(argv[1])==8){//test CurveDerivsAlg2
        iMatrix matrix = ReadDataFile_Matrix("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/Ej1.txt");
        //cout<<"abierto"<<endl;
        matrix.PrintMatrix();
        std::vector<double> Nurbs={0,0,0,0,0.25,0.5,0.75,1,1,1,1};
        int p=3;
        int dim=2;
        double N_Steps = 100000;

        std::vector<iMatrix<double>> sol(N_Steps+1);
        iMatrix<double> eval(matrix.GetNumRows(),std::min(dim,p)+1);
        double Step=1/(N_Steps);
        auto start = std::chrono::high_resolution_clock::now();
        for(unsigned int i = 0; i<N_Steps+1; i++){
            eval=CurveDerivsAlg2(Nurbs.size()-2-p,p,Nurbs,matrix, i*Step,dim-1);
            //eval=CurveDerivsAlg1(Nurbs.size()-2-p,p,Nurbs,matrix,i*Step);
            sol[i]=eval;
            //sol[i].PrintMatrix();
            //eval.PrintMatrix();
            //std::cout<< "----------------" <<std::endl;
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    
        std::cout << "Tiempo tomado: " << duration.count() << " nanosegundos\n";
        iMatrix<double> Derivadaiesima(2, N_Steps+1);
        for(unsigned int iter=0; iter<N_Steps+1;iter++){
            Derivadaiesima.SetCol(iter,sol[iter].GetRow(1));
        }
        Write_Curve(Derivadaiesima,"C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/Curvas/CurvaEj1DersAlg2.txt");
        std::cout<<"fin 8"<<endl;
    }

    else if(stoi(argv[1])==9){//test SurfacePoint
        std::vector<iMatrix<double>> CtrlPts = ReadDataFile_CtrlPts("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/Ej2.txt");
        for(unsigned int i=0; i<CtrlPts.size();i++){
            CtrlPts[i].PrintMatrix();
            std::cout<<"-------------------" << std::endl;
        }
        std::vector<double> V1={0,0,0,0.5,0.5,1,1,1};
        std::vector<double> V2={0,0,0,0,1,1,1,1};
        int p1=2;
        int p2=3;
        double N_Steps_1 = 100;
        double N_Steps_2 = 100;
        //double max_1=4;
        //double max_2=5;
        double Step_1=1/N_Steps_1;
        double Step_2=1/N_Steps_2;

        std::vector<iMatrix<double>> eval(N_Steps_1+1);
        for(unsigned int i=0; i<N_Steps_1+1;i++){
            
            iMatrix<double> auxMatrix(3,N_Steps_2+1);
            for(unsigned int j=0; j< N_Steps_2+1;j++){
                
                std::vector<double> aux = SurfacePoint(V1.size()-2-p1,p1,V1, V2.size()-2-p2,p2,V2,CtrlPts,i*Step_1,j*Step_2);
                for(unsigned int k=0; k<3;k++){
                    auxMatrix(k,j)=aux[k];
            }
            }
            eval[i]=auxMatrix;
        }

        Write_BezierSurface(eval,"C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/Superficies/SuperficieEj2" );
        Write_CtrlPts(CtrlPts,"C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/Ej2");

        std::cout << "fin 9";
    }

    if(stoi(argv[1])==10){ //Test CurvePointRational
        iMatrix CtrlPts = ReadDataFile_Matrix("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/Ej1.txt");        
        CtrlPts.PrintMatrix();
        std::vector<double> KnotVector={0,0,0,0.2,0.4,0.6,0.8,1,1,1};
        std::vector<double> Weights={1,1,0.5,1,1,1,1};
        int p=2;
        int dim=2;
        double N_Steps = 1000;
        iMatrix<double> sol(dim,N_Steps+1);
        std::vector<double> eval(dim);
        double Step=1/(N_Steps);
        std::cout<< "--------------------------------------------" << std::endl;
        iMatrix WeightedCtrlPts=WeightCtrlPts(CtrlPts,Weights);
        WeightedCtrlPts.PrintMatrix();
        for(unsigned int i = 0; i<N_Steps+1; i++){
            //std::cout << i << std::endl;
            eval=CurvePointRational(KnotVector.size()-2-p,p,KnotVector,WeightedCtrlPts,i*Step);
            //std::cout << eval[0] << "|" << eval[1] << std::endl;
            sol.SetCol(i,eval);
        }
        Write_Curve(sol,"C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/Curvas/CurvaEj1Nurbs.txt");
        std::cout<<"fin 10"<<endl;
    }
    else if(stoi(argv[1])==11){ //test RationalBasisFuns

        int p=4;
        double N_Steps = 1000;
        std::vector<double> KnotVector={0.000, 0.000, 0.000, 0.000, 0.000, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350, 0.400, 0.450, 0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.800, 0.850, 0.900, 0.950, 1.000, 1.000, 1.000, 1.000, 1.000};
        std::vector<double> WeightVector={1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000};
        iMatrix<double> sol(KnotVector.size()-p-1,N_Steps+1);
        //std::vector<double> Nurb(NurbsVector.size());
        double Step=1/(N_Steps);
        for(unsigned int i = 0; i<N_Steps+1; i++){
            double t=Step*i;
            //std::cout << "Iteración iniciada " << i<< std::endl;
            int span = FindSpan(KnotVector, t, p, KnotVector.size()-2-p);
            std::vector<double> eval= RationalBasisFuns(span, t , p, KnotVector, WeightVector);
            for(int j = span-p; j<=span;j++){
                if(j>=0) {
                    sol(j,i)=eval[j-span+p];                
                }
                else{
                    std::cout<<"Función base con indices " << j<<" "<< p << " no ha sido dibujada en el punto t="<< t << std::endl;
                    std::cout<<"Valor de la función previa: "<< eval[j-span+p]<< std::endl;
                }
                //std::cout<< "j: " << j<< " span " << span << " t:" << t << " i: " << i<<std::endl;
                
            }
            //std::cout<< " span " << span << " t:" << t <<std::endl;
            //std::cout << "Iteración completada " << i<< std::endl;
        }
        Write_BasisFunctions(sol, "C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/RationalBasisFuncts");





        std::cout<<"fin 11"<<endl;


    }


    else if(stoi(argv[1])==12){ //test RationalBasisFuns

        int p=2;
        double N_Steps = 100;
        std::vector<double>KnotVector= {0,0,0,1,2,3,3,3};
        std::vector<double>WeightVector= {1,4,1,1,1};
        iMatrix<double> sol(KnotVector.size()-p-1,N_Steps+1);
        //std::vector<double> Nurb(NurbsVector.size());
        double Step=1/(N_Steps);
        for(unsigned int i = 0; i<N_Steps+1; i++){
            double t=Step*i*3;
            //std::cout << "Iteración iniciada " << i<< std::endl;
            int span = FindSpan(KnotVector, t, p, KnotVector.size()-2-p);
            std::vector<double> eval= BasisFuns(span, t , p, KnotVector);
            std::cout<< "j sera " << span-p<< std::endl;
            for(int j = span-p; j<=span;j++){
                if(j>=0) {
                    sol(j,i)=eval[j-span+p];                    
                    std::cout<<WeightVector[j]<<"|";
                }
                else{
                    std::cout<<"Función base con indices " << j<<" "<< p << " no ha sido dibujada en el punto t="<< t << std::endl;
                    std::cout<<"Valor de la función previa: "<< eval[j-span+p]<< std::endl;
                }
                //std::cout<< "j: " << j<< " span " << span << " t:" << t << " i: " << i<<std::endl;
                
            }
            std::cout<< " span " << span << " t:" << t <<std::endl;
            //std::cout << "Iteración completada " << i<< std::endl;
        }
        Write_BasisFunctions(sol, "C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/BasisFuncts2");





        std::cout<<"fin 12"<<endl;


    }


     if(stoi(argv[1])==13){ //Test CurvePointRational
        iMatrix CtrlPts = ReadDataFile_Matrix("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/EjPruebaFEM.txt");
        //cout<<"abierto"<<endl;
        CtrlPts.PrintMatrix();
                std::cout<< "--------------------------------------------" << std::endl;
        double t = 0.885;
        std::vector<double> KnotVector={0,0,0,1,1,1};
        std::vector<double> Weights={1,1,1};
        iMatrix WeightedCtrlPtsPreDrop=WeightCtrlPts(CtrlPts,Weights);
        int p=2;
        int n= KnotVector.size()-p-2;
        std::vector<double> Insertions={0.2,0.4,0.6,0.8};
        iMatrix<double> NewCtrlPts=RefineKnotVectCurve(n,p,KnotVector, WeightedCtrlPtsPreDrop, Insertions, Insertions.size()-1);
        std::cout<< "----------------NewCtrlPtsW--------------" << std::endl; 
        NewCtrlPts.PrintMatrix();  
        iMatrix<double> WeightsMat(1,NewCtrlPts.GetNumCols());
        WeightsMat.SetRow(0,NewCtrlPts.GetRow(2));
        std::cout<< "----------------WeightsMat--------------" << std::endl; 
        WeightsMat.PrintMatrix();
    
        
        int dim=2;
        double N_Steps = 100;
        iMatrix<double> sol(dim,N_Steps+1);
        std::vector<double> eval(dim);
        double Step=1/(N_Steps);
        std::cout<< "--------------WeightedCtrlPts-----------------"<< std::endl;
        
        iMatrix<double> WeightedCtrlPts(NewCtrlPts.GetNumRows()-1,NewCtrlPts.GetNumCols());
        for(unsigned int i=0; i< WeightedCtrlPts.GetNumRows();i++){
            WeightedCtrlPts.SetRow(i,NewCtrlPts.GetRow(i));
        }
        WeightedCtrlPts.PrintMatrix();
        iMatrix<double> Allders= CurveDerivsAlg1(KnotVector.size()-p-2, p, KnotVector, NewCtrlPts, t, 1);
        iMatrix<double> Aders = CurveDerivsAlg1(KnotVector.size()-p-2, p, KnotVector, WeightedCtrlPts, t, 1);
        //std::cout << "p " << p<< " d " << 2 << std::endl;
        iMatrix<double> wders = CurveDerivsAlg1(KnotVector.size()-p-2, p, KnotVector, WeightsMat, t, 1);
        //std::cout << "p " << p<< " d " << 2 << std::endl;
        std::cout<<"----------Allders----------"<<std::endl;
        Allders.PrintMatrix();
        std::cout<<"----------Aders----------"<<std::endl;
        Aders.PrintMatrix();
        std::cout<<"----------wders-------------"<<std::endl;
        wders.PrintMatrix();
        std::cout<<"-----------CW---------------"<<std::endl;


        iMatrix<double> CW= RatCurveDerivs(Aders, wders, 1);
        CW.PrintMatrix();


        iMatrix<double> Aders2=Allders.GetSubMat(0,1,0,1);

        iMatrix<double> wders2=Allders.GetSubMat(0,1,2,2);

        std::cout<<"----------Aders2----------"<<std::endl;
        Aders2.PrintMatrix();
        std::cout<<"----------wders2-------------"<<std::endl;
        wders2.PrintMatrix();

        


        std::cout<<"fin 13"<<endl;
    }

    if(stoi(argv[1])==14){ //Test SurfaceRational
        
        std::vector<iMatrix<double>> CtrlPts = ReadDataFile_CtrlPts("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/Ej2.txt");
        std::cout<<"----CtrlPts----" << std::endl;
        for(unsigned int i=0; i<CtrlPts.size();i++){
            CtrlPts[i].PrintMatrix();
            std::cout<<"-------------------" << std::endl;
        }
        iMatrix<double> Weights =ReadDataFile_Matrix("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Weights/Ej2Weights.txt");
        std::cout<<"---Matriz_Pesos---" << std::endl;
        Weights.PrintMatrix();
        std::cout<<"-----CtrlPtsW-----" << std::endl;
        std::vector<iMatrix<double>> CtrlPtsW=WeightCtrlPtsSurf(CtrlPts,Weights);
        for(unsigned int i=0; i<CtrlPtsW.size();i++){
            CtrlPtsW[i].PrintMatrix();
            std::cout<<"-------------------" << std::endl;
        }

        std::vector<double> V1={0,0,0,0.5,0.5,1,1,1};
        std::vector<double> V2={0,0,0,0,1,1,1,1};
        int p1=2;
        int p2=3;
        double N_Steps_1 = 100;
        double N_Steps_2 = 100;

        double Step_1=1/N_Steps_1;
        double Step_2=1/N_Steps_2;

        std::vector<iMatrix<double>> eval(N_Steps_1+1);
        for(unsigned int i=0; i<N_Steps_1+1;i++){
            
            iMatrix<double> auxMatrix(3,N_Steps_2+1);
            for(unsigned int j=0; j< N_Steps_2+1;j++){
                
                std::vector<double> aux = SurfacePointRational(V1.size()-2-p1,p1,V1, V2.size()-2-p2,p2,V2,CtrlPtsW,i*Step_1,j*Step_2);
                for(unsigned int k=0; k<3;k++){
                    auxMatrix(k,j)=aux[k];
            }
            }
            eval[i]=auxMatrix;
        }

        Write_BezierSurface(eval,"C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/Superficies/SuperficieRationalEj2" );
        Write_CtrlPts(CtrlPts,"C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/Ej2");

    


        std::cout<<"fin 14"<<endl;
    }
    if(stoi(argv[1])==15){ //Test RatSurfaceDers
        
        double t1= 0.8f;
        double t2=0.8f;
        int p1=2;
        int p2=3;
        std::vector<double> V1={0,0,0,0.5,0.5,1,1,1};
        std::vector<double> V2={0,0,0,0,1,1,1,1};
        int MaxD=2;
        std::vector<iMatrix<double>> CtrlPts = ReadDataFile_CtrlPts("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/Ej2.txt");
        std::cout<<"----CtrlPts----" << std::endl;
        for(unsigned int i=0; i<CtrlPts.size();i++){
            CtrlPts[i].PrintMatrix();
            std::cout<<"-------------------" << std::endl;
        }
        iMatrix<double> Weights =ReadDataFile_Matrix("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Weights/Ej2Weights.txt");
        std::cout<<"---Matriz_Pesos---" << std::endl;
        Weights.PrintMatrix();
        std::cout<<"-----CtrlPtsW-----" << std::endl;
        std::vector<iMatrix<double>> CtrlPtsW=WeightCtrlPtsSurf(CtrlPts,Weights);
        
        for(unsigned int i=0; i<CtrlPtsW.size();i++){
            CtrlPtsW[i].PrintMatrix();
            std::cout<<"-------------------" << std::endl;
        }  
        
        std::vector<iMatrix<double>> CtrlPtsW2(CtrlPtsW.size());
        int dim=CtrlPtsW[0].GetNumRows()-1;
        for(unsigned int i=0; i<CtrlPtsW.size();i++){
            iMatrix<double> aux(dim, CtrlPtsW[0].GetNumCols());
            for(unsigned int j=0; j<dim;j++){
                aux.SetRow(j,CtrlPtsW[i].GetRow(j));
                
            }
            CtrlPtsW2[i]=aux;
        }



        iMatrix<std::vector<double>> Aders = SurfaceDerivsAlg1(V1.size()-2-p1,p1,V1,V2.size()-2-p2,p2,V2,CtrlPtsW2, t1, t2, MaxD);
        std::cout<<"--------Aders------" << std::endl;
        for(unsigned int i=0; i<=MaxD; i++){
            for(unsigned int j=0; j<=MaxD; j++){
                for(unsigned int k=0; k<Aders(0,0).size();k++){
                    std::cout << Aders(i,j)[k]<< " ";
                }
                std::cout<<"|";
            }
            std::cout << std::endl;
        }


        std::cout<<"------Weights2-----" << std::endl;
        
        std::vector<iMatrix<double>> Weights2(Weights.GetNumRows());
        dim=1;
        for(unsigned int i=0; i<CtrlPtsW.size();i++){
            iMatrix<double> aux(1, CtrlPtsW[0].GetNumCols());

                aux.SetRow(0,CtrlPtsW[i].GetRow(CtrlPtsW[0].GetNumRows()-1));
                
            Weights2[i]=aux;
            Weights2[i].PrintMatrix();
                    std::cout<<"-----------" << std::endl;
        }
        std::cout<<"------Wders-----" << std::endl;
        iMatrix<std::vector<double>> Wders = SurfaceDerivsAlg1(V1.size()-2-p1,p1,V1,V2.size()-2-p2,p2,V2,Weights2, t1, t2, MaxD);
        iMatrix<double> WdersAux(MaxD+1,MaxD+1);
        for(unsigned int i=0; i<=MaxD; i++){
            for(unsigned int j=0; j<=MaxD; j++){
                //for(unsigned int k=0; k<Wders(0,0).size();k++){
                    std::cout << Wders(i,j)[0]<< " ";
                    WdersAux(i,j)=Wders(i,j)[0];
               //}
                std::cout<<"|";
            }
            std::cout << std::endl;
        }
        WdersAux.PrintMatrix();

        iMatrix<std::vector<double>> SurfRatDers= RatSurfaceDerivs(Aders, WdersAux, MaxD);
        std::cout<< "numRows: " <<SurfRatDers.GetNumRows()<< ", numCols: "<<SurfRatDers.GetNumCols()<< std::endl;
        for(unsigned int i=0; i<= MaxD; i++){
            for(unsigned int j=0; j<= MaxD; j++){
                for(unsigned int k=0; k<SurfRatDers(0,0).size();k++){
                    if(i+j<=2){
                         std::cout<<SurfRatDers(i,j)[k]<<" ";
                    }
                    else{
                        std::cout<<"0 ";
                    }
                   
                }
                std::cout<<"|";
            }
            std::cout<< std::endl;
        }




        std::cout<<"fin 15"<<endl;    
    }

    if(stoi(argv[1])==16){ //Test Geometric Algorithms 1
        
        iMatrix CtrlPts = ReadDataFile_Matrix("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/Ej1.txt");
        //cout<<"abierto"<<endl;
        CtrlPts.PrintMatrix();
                std::cout<< "--------------------------------------------" << std::endl;
        std::vector<double> KnotVector={0,0,0,0.2,0.4,0.6,0.8,1,1,1};
        double element = 0.2f;
        int index=findInsertionIndex(KnotVector,element);
        int r = 1;
        std::vector<double> Weights={1,1,0.5,1,1,1,1};
        iMatrix<double> WeightsMat(1,Weights.size());
        WeightsMat.SetRow(0,Weights);
        WeightsMat.PrintMatrix();
        std::cout<< "--------------------------------------------"<< std::endl;
        iMatrix WeightedCtrlPts=WeightCtrlPts(CtrlPts,Weights);
        WeightedCtrlPts.PrintMatrix();
        std::cout<< "--------------------------------------------"<< std::endl;
        int p=2;
        int dim=2;
        iMatrix<double> NewCtrlPtsW= CurveKnotsIns(KnotVector.size()-2-p, p, KnotVector,WeightedCtrlPts, element, index, 0, r);
        NewCtrlPtsW.PrintMatrix();
        double N_Steps = 100;
        iMatrix<double> sol(dim,N_Steps+1);
        std::vector<double> eval(dim);
        double Step=1/(N_Steps);
        std::cout<< "--------------------------------------------" << std::endl;
        for(unsigned int i=0; i<KnotVector.size();i++){
            std::cout << KnotVector[i] << ", ";
        }
        std::cout<< std::endl<< "--------------------------------------------"<< std::endl;
        for(unsigned int i = 0; i<N_Steps+1; i++){
            //std::cout << i << std::endl;
            eval=CurvePointRational(KnotVector.size()-2-p,p,KnotVector,NewCtrlPtsW,i*Step);
            //std::cout << eval[0] << "|" << eval[1] << std::endl;
            sol.SetCol(i,eval);
        }
        Write_Curve(sol,"C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/Curvas/CurvaEj1Nurbs.txt");
        std::cout<< "fin 16";
    }
    if(stoi(argv[1])==17){ //Test CurvePntByCornerCut
        
        iMatrix CtrlPts = ReadDataFile_Matrix("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/Ej1.txt");
        std::vector<double> KnotVector={0,0,0,0.2,0.4,0.6,0.8,1,1,1};
        
        

        std::vector<double> Weights={1,1,0.5,1,1,1,1};
        iMatrix<double> WeightsMat(1,Weights.size());
        WeightsMat.SetRow(0,Weights);
        iMatrix WeightedCtrlPts=WeightCtrlPts(CtrlPts,Weights);
        int p=2;
        int dim=2;
        /*
        int span1 = FindSpan(KnotVector, element, p, KnotVector.size()-2-p);
        int mult1= countOccurrences(KnotVector,element);
        int span2, mult2;
        FindSpanMult(KnotVector, element, p, KnotVector.size()-2-p, span2, mult2);
        std::cout << "span1: "<< span1 << ", mult1: " << mult1<< std::endl;
        std::cout << "span2: "<< span2 << ", mult2: " << mult2<< std::endl; 
        */
        double N_Steps = 100;
        iMatrix<double> sol(dim,N_Steps+1);
        std::vector<double> eval(dim);
        double Step=1/(N_Steps);
        for(unsigned int i = 0; i<N_Steps+1; i++){
            //std::cout << i << std::endl;
            eval=CurvePntByCornerCut(KnotVector.size()-2-p,p,KnotVector,WeightedCtrlPts,i*Step);
            //std::cout << eval[0] << "|" << eval[1] << std::endl;
            sol.SetCol(i,eval);
        }
        Write_Curve(sol,"C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/Curvas/CurvaEj1CornerCut.txt");
        std::cout<<"fin 17";
    }

    if(stoi(argv[1])==18){//test SurfaceKnotIns
        std::vector<iMatrix<double>> CtrlPts = ReadDataFile_CtrlPts("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/Ej2.txt");
        std::cout<<"----CtrlPts----" << std::endl;
        for(unsigned int i=0; i<CtrlPts.size();i++){
            CtrlPts[i].PrintMatrix();
            std::cout<<"-------------------" << std::endl;
        }
        iMatrix<double> Weights =ReadDataFile_Matrix("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Weights/Ej2Weights.txt");
        std::cout<<"---Matriz_Pesos---" << std::endl;
        Weights.PrintMatrix();
        std::cout<<"-----CtrlPtsW-----" << std::endl;
        std::vector<iMatrix<double>> CtrlPtsW=WeightCtrlPtsSurf(CtrlPts,Weights);
        for(unsigned int i=0; i<CtrlPtsW.size();i++){
            CtrlPtsW[i].PrintMatrix();
            std::cout<<"-------------------" << std::endl;
        }

        std::vector<double> V1={0,0,0,0.5,0.5,1,1,1};
        std::vector<double> V2={0,0,0,0,1,1,1,1};
        int p1=2;
        int p2=3;

        std::vector<iMatrix<double>> NewCtrlPtsW= SurfaceKnotIns(V1.size()-p1-2, p1, V1, V2.size()-p2-2,p2, V2, CtrlPtsW, false, 0.2f, findInsertionIndex(V2,0.2f) ,1,0);
        std::cout<<"------NewCtrlPtsW-----" << std::endl;
        for(unsigned int i=0; i<NewCtrlPtsW.size();i++){
            NewCtrlPtsW[i].PrintMatrix();
            std::cout<<"-------------------" << std::endl;
        }
        std::cout<<"------KnotVector1-----" << std::endl;
        for(unsigned int i=0; i<V1.size();i++){
            
            std::cout<<V1[i]<< ", ";
        }
        std::cout<<std::endl;
        std::cout<<"------KnotVector2-----" << std::endl;
        for(unsigned int i=0; i<V2.size();i++){
            
            std::cout<<V2[i]<< ", ";
        }
        double N_Steps_1 = 100;
        double N_Steps_2 = 100;

        double Step_1=1/N_Steps_1;
        double Step_2=1/N_Steps_2;

        std::vector<iMatrix<double>> eval(N_Steps_1+1);
        for(unsigned int i=0; i<N_Steps_1+1;i++){
            
            iMatrix<double> auxMatrix(3,N_Steps_2+1);
            for(unsigned int j=0; j< N_Steps_2+1;j++){
                
                std::vector<double> aux = SurfacePointRational(V1.size()-2-p1,p1,V1, V2.size()-2-p2,p2,V2,NewCtrlPtsW,i*Step_1,j*Step_2);
                for(unsigned int k=0; k<3;k++){
                    auxMatrix(k,j)=aux[k];
            }
            }
            eval[i]=auxMatrix;
        }
        std::cout<< std::endl;
        Write_BezierSurface(eval,"C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/Superficies/SuperficieRationalEj2Insertion" );
        std::vector<iMatrix<double>> NewCtrlPts = UnWeightCtrlPtsSurf(NewCtrlPtsW);
        std::cout<<"------NewCtrlPts-----" << std::endl;
        for(unsigned int i=0; i<NewCtrlPts.size();i++){
            NewCtrlPts[i].PrintMatrix();
            std::cout<<"-------------------" << std::endl;
        }
        Write_CtrlPts(NewCtrlPts,"C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/Ej2Insertion");
        std::cout<<"fin 18";
    }

    if(stoi(argv[1])==19){//test RefineKnotVectCurve
        iMatrix CtrlPts = ReadDataFile_Matrix("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/Ej1.txt"); 
        std::cout<< "----------------CtrlPts--------------" << std::endl;       
        CtrlPts.PrintMatrix();
        std::vector<double> KnotVector={0,0,0,0.2,0.4,0.6,0.8,1,1,1};
        std::vector<double> Insertions={0.1,0.3,0.5,0.7,0.9};
        std::vector<double> Weights={1,1,0.5,1,1,1,1};
        int p=2;
        int dim=2;
        double N_Steps = 100;
        iMatrix<double> sol(dim,N_Steps+1);
        std::vector<double> eval(dim);
        double Step=1/(N_Steps);
        std::cout<< "----------------WeightedCtrlPts--------------" << std::endl;
        iMatrix WeightedCtrlPts=WeightCtrlPts(CtrlPts,Weights);
        WeightedCtrlPts.PrintMatrix();
        iMatrix<double> NewCtrlPtsW= RefineKnotVectCurve(KnotVector.size()-p-2,p,KnotVector,WeightedCtrlPts, Insertions, Insertions.size()-1);
        std::cout<< "----------------NewCtrlPtsW--------------" << std::endl;
        NewCtrlPtsW.PrintMatrix();
        std::cout<< "----------------NewKnotVector--------------" << std::endl;
        for(unsigned int i=0; i<KnotVector.size();i++){
            std::cout << KnotVector[i] << ", ";
        }
        std::cout<< std::endl;

        for(unsigned int i = 0; i<N_Steps+1; i++){
            //std::cout << i << std::endl;
            eval=CurvePointRational(KnotVector.size()-2-p,p,KnotVector,NewCtrlPtsW,i*Step);
            //std::cout << eval[0] << "|" << eval[1] << std::endl;
            sol.SetCol(i,eval);
        }
        Write_Curve(sol,"C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/Curvas/CurvaEj1RefineKnotVect.txt");
        std::cout<<"fin 19"<<endl;
    }
    if(stoi(argv[1])==20){//test RefineKnotVectSurf
         std::vector<iMatrix<double>> CtrlPts = ReadDataFile_CtrlPts("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/Ej2.txt");
        std::cout<<"----CtrlPts----" << std::endl;
        for(unsigned int i=0; i<CtrlPts.size();i++){
            CtrlPts[i].PrintMatrix();
            std::cout<<"-------------------" << std::endl;
        }
        iMatrix<double> Weights =ReadDataFile_Matrix("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Weights/Ej2Weights.txt");
        std::cout<<"---Matriz_Pesos---" << std::endl;
        Weights.PrintMatrix();
        std::cout<<"-----CtrlPtsW-----" << std::endl;
        std::vector<iMatrix<double>> CtrlPtsW=WeightCtrlPtsSurf(CtrlPts,Weights);
        for(unsigned int i=0; i<CtrlPtsW.size();i++){
            CtrlPtsW[i].PrintMatrix();
            std::cout<<"-------------------" << std::endl;
        }

        std::vector<double> Insertions={0.2,0.4,0.6,0.8};
        std::vector<double> V1={0,0,0,0.5,0.5,1,1,1};
        std::vector<double> V2={0,0,0,0,1,1,1,1};
        int p1=2;
        int p2=3;

        std::vector<iMatrix<double>> NewCtrlPtsW= RefineKnotVectSurface(V1.size()-p1-2, p1, V1, V2.size()-p2-2,p2, V2, CtrlPtsW, false, Insertions, Insertions.size()-1);
        NewCtrlPtsW= RefineKnotVectSurface(V1.size()-p1-2, p1, V1, V2.size()-p2-2,p2, V2, NewCtrlPtsW, true, Insertions, Insertions.size()-1);
        std::cout<<"------NewCtrlPtsW-----" << std::endl;
        for(unsigned int i=0; i<NewCtrlPtsW.size();i++){
            NewCtrlPtsW[i].PrintMatrix();
            std::cout<<"-------------------" << std::endl;
        }
        std::cout<<"------KnotVector1-----" << std::endl;
        for(unsigned int i=0; i<V1.size();i++){
            
            std::cout<<V1[i]<< ", ";
        }
        std::cout<<std::endl;
        std::cout<<"------KnotVector2-----" << std::endl;
        for(unsigned int i=0; i<V2.size();i++){
            
            std::cout<<V2[i]<< ", ";
        }
        std::cout<< std::endl;
        double N_Steps_1 = 100;
        double N_Steps_2 = 100;

        double Step_1=1/N_Steps_1;
        double Step_2=1/N_Steps_2;

        std::vector<iMatrix<double>> eval(N_Steps_1+1);
        for(unsigned int i=0; i<N_Steps_1+1;i++){
            
            iMatrix<double> auxMatrix(3,N_Steps_2+1);
            for(unsigned int j=0; j< N_Steps_2+1;j++){
                
                std::vector<double> aux = SurfacePointRational(V1.size()-2-p1,p1,V1, V2.size()-2-p2,p2,V2,NewCtrlPtsW,i*Step_1,j*Step_2);
                for(unsigned int k=0; k<3;k++){
                    auxMatrix(k,j)=aux[k];
            }
            }
            eval[i]=auxMatrix;
        }
        std::cout<< std::endl;
        Write_BezierSurface(eval,"C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/Superficies/SuperficieRationalEj2RefineKnotVectSurf" );
        std::vector<iMatrix<double>> NewCtrlPts = UnWeightCtrlPtsSurf(NewCtrlPtsW);
        std::cout<<"------NewCtrlPts-----" << std::endl;
        for(unsigned int i=0; i<NewCtrlPts.size();i++){
            NewCtrlPts[i].PrintMatrix();
            std::cout<<"-------------------" << std::endl;
        }
        Write_CtrlPts(NewCtrlPts,"C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/SuperficieRationalEj2RefineKnotVectSurf");
        std::cout<<"fin 20";
    }

    if(stoi(argv[1])==21){ //Test SurfaceRational cilindro
        
        std::vector<iMatrix<double>> CtrlPts = ReadDataFile_CtrlPts("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/Ej3.txt");
        std::cout<<"----CtrlPts----" << std::endl;
        for(unsigned int i=0; i<CtrlPts.size();i++){
            CtrlPts[i].PrintMatrix();
            std::cout<<"-------------------" << std::endl;
        }
        iMatrix<double> Weights =ReadDataFile_Matrix("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Weights/Ej3Weights.txt");
        std::cout<<"---Matriz_Pesos---" << std::endl;
        Weights.PrintMatrix();
        std::cout<<"-----CtrlPtsW-----" << std::endl;
        std::vector<iMatrix<double>> CtrlPtsW=WeightCtrlPtsSurf(CtrlPts,Weights);
        for(unsigned int i=0; i<CtrlPtsW.size();i++){
            CtrlPtsW[i].PrintMatrix();
            std::cout<<"-------------------" << std::endl;
        }

        std::vector<double> V1={0,0,0,0.25,0.25,0.5,0.5,0.75,0.75,1,1,1};
        std::vector<double> V2={0,0,0.250,0.500,0.750,1,1};
        int p1=2;
        int p2=1;
        double N_Steps_1 = 100;
        double N_Steps_2 = 100;

        double Step_1=1/N_Steps_1;
        double Step_2=1/N_Steps_2;

        std::vector<iMatrix<double>> eval(N_Steps_1+1);
        for(unsigned int i=0; i<N_Steps_1+1;i++){
            
            iMatrix<double> auxMatrix(3,N_Steps_2+1);
            for(unsigned int j=0; j< N_Steps_2+1;j++){
                
                std::vector<double> aux = SurfacePointRational(V1.size()-2-p1,p1,V1, V2.size()-2-p2,p2,V2,CtrlPtsW,i*Step_1,j*Step_2);
                for(unsigned int k=0; k<3;k++){
                    auxMatrix(k,j)=aux[k];
            }
            }
            eval[i]=auxMatrix;
        }

        Write_BezierSurface(eval,"C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/Superficies/SuperficieRationalEj3" );
        Write_CtrlPts(CtrlPts,"C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/Ej3");

        std::vector<double> Insertions1={0.1,0.2,0.3,0.4,0.6,0.7,0.8,0.9};
        std::vector<iMatrix<double>> NewCtrlPtsW= RefineKnotVectSurface(V1.size()-p1-2, p1, V1, V2.size()-p2-2,p2, V2, CtrlPtsW, true, Insertions1, Insertions1.size()-1);
        std::vector<double> Insertions2={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
        NewCtrlPtsW= RefineKnotVectSurface(V1.size()-p1-2, p1, V1, V2.size()-p2-2,p2, V2, NewCtrlPtsW, false, Insertions2, Insertions2.size()-1);
        std::cout<<"------NewCtrlPtsW-----" << std::endl;
        for(unsigned int i=0; i<NewCtrlPtsW.size();i++){
            NewCtrlPtsW[i].PrintMatrix();
            std::cout<<"-------------------" << std::endl;
        }
        std::cout<<"------KnotVector1-----" << std::endl;
        for(unsigned int i=0; i<V1.size();i++){
            
            std::cout<<V1[i]<< ", ";
        }
        std::cout<<std::endl;
        std::cout<<"------KnotVector2-----" << std::endl;
        for(unsigned int i=0; i<V2.size();i++){
            
            std::cout<<V2[i]<< ", ";
        }
         std::cout<< std::endl;
        for(unsigned int i=0; i<N_Steps_1+1;i++){
            
            iMatrix<double> auxMatrix(3,N_Steps_2+1);
            for(unsigned int j=0; j< N_Steps_2+1;j++){
                
                std::vector<double> aux = SurfacePointRational(V1.size()-2-p1,p1,V1, V2.size()-2-p2,p2,V2,NewCtrlPtsW,i*Step_1,j*Step_2);
                for(unsigned int k=0; k<3;k++){
                    auxMatrix(k,j)=aux[k];
            }
            }
            eval[i]=auxMatrix;
        }
    
        std::vector<iMatrix<double>> NewCtrlPts = UnWeightCtrlPtsSurf(NewCtrlPtsW);
        std::cout<<"------NewCtrlPts-----" << std::endl;
        for(unsigned int i=0; i<NewCtrlPts.size();i++){
            NewCtrlPts[i].PrintMatrix();
            std::cout<<"-------------------" << std::endl;
        }
        
        Write_BezierSurface(eval,"C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/Superficies/SuperficieRationalEj3Refined" );
        Write_CtrlPts(NewCtrlPts,"C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/Ej3Refined");


        std::cout<<"fin 21"<<endl;
    }
    /*
    if(stoi(argv[1])==22){

        SparseTensor<double> tensor(3); 

        tensor({0,)=0, 0}, 0.5);
        tensor({0,)=1, 0}, 3.5);
        tensor({0,)=1, 1}, 3.4);
        tensor({0,)=2, 1}, 7.7);
        tensor({2,)=3, 2}, -1.2);
        tensor({2,)=1, 2}, -2.4);
        tensor({0,)=1, 3}, 0.0); // limina el valor
        tensor({3,)=3, 3}, 322.5);
        tensor({2,)=1, 3}, 99);

    //std::cout << "Valor en (2, 3): " << tensor.get({2, 3}) << std::endl;

        tensor.forEach([](const std::vector<int>& idx, double val) {
            std::cout << "(";
            for (size_t i = 0; i < idx.size(); ++i) {
                std::cout << idx[i] << (i + 1 < idx.size() ? ", " : "");
            }
            std::cout << ") = " << val << std::endl;
        });
        
        tensor.forEachComp(0, 2, [](const std::vector<int>& idx, double val) {
            std::cout << "(";
            for (size_t i = 0; i < idx.size(); ++i) {
                std::cout << idx[i] << (i + 1 < idx.size() ? ", " : "");
            }
            std::cout << "), val = " << val << std::endl;
        });

        Eigen::SparseMatrix<double> mat = tensor.toEigenSparse();

        std::cout<<"Numero de elemntos no nulos= " << mat.nonZeros()<<std::endl;
        std::cout<<"Numero de filas= " << mat.rows()<<std::endl;
        std::cout<<"Numero de columnas= " << mat.cols()<<std::endl;

        for(unsigned int i=0; i<mat.rows(); i++){
            for(unsigned int j=0; j<mat.cols();j++){
                std::cout<< mat.coeff(i,j) << ", ";
            }
            std::cout << std::endl;
        }

        std::cout<<"fin 22"<< std::endl;

    }
    */
    if(stoi(argv[1])==23){//test Legendre plynomial
        int n = 5; 
        Eigen::VectorXd nodes, weights;

        legendre_pol(n, nodes, weights);

        std::cout << "Nodos (raices):\n" << nodes.transpose() << "\n";
        std::cout << "Pesos:\n" << weights.transpose() << "\n"; 

        std::cout<<"fin 23"<< std::endl;
        
    }

    if(stoi(argv[1])==24){//test RatDersBasisFuns
        int p=3;
        std::vector<double> Nurbs={0,0,0,0,2,4,4,6,6,6,8,8,8,8};
        std::vector<double> weights={1,1,2,0.5,0.5,7,1,1,1};


        double N_Steps = 100;
        std::vector<double> eval(Nurbs.size()-p);
        iMatrix<double> ders(p+1,p+1);
        iMatrix<double> sol(Nurbs.size()-p,N_Steps+1);
        //std::vector<double> Nurb(NurbsVector.size());
        double Step=1/(N_Steps);
        auto start = std::chrono::high_resolution_clock::now();
        for(unsigned int i = 0; i<N_Steps+1; i++){
            
            double t=Step*i*8;
            //std::cout << "Iteración iniciada " << i<< std::endl;
            int span = FindSpan(Nurbs, t, p, Nurbs.size()-2-p);
            ders = RatDersBasisFuns(span,t,p,1,Nurbs, weights);
            ders.PrintMatrix();
            std::cout<<"---------------------"<<std::endl;
            for(int j = span-p; j<=span;j++){
                if(j>=0) {
                    sol(j,i)=ders(1,j-span+p);
                }
                else{
                    std::cout<<"Función base con indices " << j<<" "<< p << " no ha sido dibujada en el punto t="<< t << std::endl;
                    std::cout<<"Valor de la función previa: "<< eval[j-span+p]<< std::endl;
                }
                //std::cout<< "j: " << j<< " span " << span << " t:" << t << " i: " << i<<std::endl;
                
            }
        }
        
        Write_BasisFunctions(sol, "C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/RatBasisFuncts_Derivates");
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    
        std::cout << "Tiempo tomado: " << duration.count() << " nanosegundos\n";
    
        std::cout<<"fin 24"<<endl;
    }

    if(stoi(argv[1])==25){//test 1D problems (Unoptimized legacy functions)
        /*
        bool PrintMatixes=false;
        bool PrintIntermediateValues=false;

        
        int p=4;
        double nElements = 40;
        //variables que no cambian en cada ejecución del bucle o que se declaran fuera para ir actualizandolas
        std::cout<< "----------------CtrlPts--------------" << std::endl;
        iMatrix<double> CtrlPts = ReadDataFile_Matrix("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/EjPruebaFEM.txt");        
        CtrlPts.PrintMatrix();
        std::vector<double> KnotVector={0,0,0,0,0,1,1,1,1,1};
        double lower_limit=KnotVector[0];
        double upper_limit=KnotVector[KnotVector.size()-1];
        std::vector<double> Weights={1,1,1,1,1};
        std::cout<< "----------WeightedCtrlPts--------------" << std::endl;
        iMatrix WeightedCtrlPts=WeightCtrlPts(CtrlPts,Weights);
        WeightedCtrlPts.PrintMatrix();
        
        int n= KnotVector.size()-p-2;
        std::cout<< "---------------KnotVectorOld-----------------" << std::endl;
        for(unsigned int i=0; i<KnotVector.size(); i++){
            std::cout << KnotVector[i] << ", ";
        }
        std::cout << std::endl;
        std::cout << "n value:" << n << std::endl;
        
        std::vector<double> Insertions=createSubIntervals(nElements);
        Insertions=removeMatches(Insertions, KnotVector);
        iMatrix<double> NewCtrlPts=RefineKnotVectCurve(n,p,KnotVector, WeightedCtrlPts, Insertions, Insertions.size()-1);
        std::cout<< "----------------NewCtrlPtsW--------------" << std::endl;   
        NewCtrlPts.PrintMatrix();
        int size = KnotVector.size()-p-1;
        auto f = [](double x) { 
            double pi=3.141592653589793;
            return pi*pi*sin(pi*x); 
        };

        auto analyticSol = [](double x){
            double pi=3.141592653589793;
            return sin(pi*x);
        };

        auto analyticSolDers = [](double x){
            double pi=3.141592653589793;
            return pi*cos(pi*x);
        };

        n = KnotVector.size()-p-2;

        std::cout<< "---------------KnotVectorNew-----------------" << std::endl;
        for(unsigned int i=0; i<KnotVector.size(); i++){
            std::cout << KnotVector[i] << ", ";
        }
        std::cout << std::endl;
        std::cout << "n value:" << n << std::endl;

        Weights=NewCtrlPts.GetRow(NewCtrlPts.GetNumRows()-1);
        std::cout<< "---------------Weights-----------------" << std::endl;
        for(unsigned int i=0; i<Weights.size(); i++){
            std::cout << Weights[i] << ", ";
        }
        
        int dim=NewCtrlPts.GetNumRows()-1;
        int numElements= NewCtrlPts.GetNumCols();
        
        int nEvals=p+2;
        
        
        Eigen::SparseMatrix<double> global(size, size);
        Eigen::VectorXd LinearForm(size);
        LinearForm.setZero(); // Inicializa en cero
        
        
        std::cout<< "---------------Intervals-----------------" << std::endl;
        std::vector<double> Intervals=removeDuplicates(KnotVector);
        for(unsigned int i=0; i<Intervals.size(); i++){
            std::cout << Intervals[i] << ", ";
        }
        std::cout << std::endl;
        Eigen::VectorXd nodes, weights;
        legendre_pol(nEvals, nodes, weights);
            
        std::cout << "Nodos (raices):\n" << nodes.transpose() << "\n";
        std::cout << "Pesos:\n" << weights.transpose() << "\n"; 
        typedef Eigen::Triplet<double> T;
        std::vector<Eigen::Triplet<double>> tripletList;

        //iteraciones 
        auto start = std::chrono::high_resolution_clock::now();
        for(unsigned int index=0; index<Intervals.size()-1; index++){
            double lowerLimit=Intervals[index], upperLimit=Intervals[index+1];
            int span = FindSpan(KnotVector, lowerLimit,p,n);
            std::vector<int> global_indices(p+1);
            for(unsigned int w=0; w<=p; w++){
                global_indices[w]=span-p+w;
            }
            if(PrintIntermediateValues){

                std::cout << "global_indices: (";
                for(unsigned int w=0; w<p; w++){
                    std::cout <<global_indices[w] << ", ";
                }
                std::cout << global_indices[p] << ")" <<std::endl;
                std::cout << "Evaluation Interval: (" << Intervals[index] << ", " << Intervals[index+1]<< ")" <<std::endl;
                std::cout<< "--------------------------------------------" << std::endl;
            }
            
            iMatrix<double> BasisFunsEvals;
            iMatrix<double> DerBasisFunsEvals;
            std::vector<double> JacobianEvals;
            std::vector<double> InverseJacobianEvals;
            std::vector<double> funcEvals;
            
            D1_element_eval(n, span,p,nEvals, lowerLimit, upperLimit, KnotVector, Weights, nodes, NewCtrlPts, 
                            f, BasisFunsEvals, DerBasisFunsEvals, JacobianEvals, InverseJacobianEvals, funcEvals);
            if(PrintIntermediateValues){
                std::cout<< "-------------BasisFunsEvals----------------" << std::endl;
                BasisFunsEvals.PrintMatrix();
                std::cout<< "-------------DerBasisFunsEvals----------------" << std::endl;
                DerBasisFunsEvals.PrintMatrix();
                std::cout<< "--------------JacobianEvals----------------" << std::endl;
                for(unsigned int i=0; i<JacobianEvals.size(); i++){
                    std::cout << JacobianEvals[i] << ", ";
                }
                std::cout << std::endl;
                std::cout<< "--------------InverseJacobianEvals----------------" << std::endl;
                for(unsigned int i=0; i<InverseJacobianEvals.size(); i++){
                    std::cout << InverseJacobianEvals[i] << ", ";
                }
                std::cout << std::endl;
                std::cout<< "--------------funcEvals----------------" << std::endl;
                for(unsigned int i=0; i<funcEvals.size(); i++){
                    std::cout << funcEvals[i] << ", ";
                }
                std::cout << std::endl;
                std::cout<< "--------------Element----------------" << std::endl;
            }

            iMatrix<double> Element=gauss_legendre_cuadrature_integral_bilinealForm(p,lowerLimit,upperLimit,nEvals, DerBasisFunsEvals, InverseJacobianEvals, weights);
            
            if(PrintIntermediateValues){
                Element.PrintMatrix();
                std::cout<< "-------------------------------------" << std::endl;
            }


            for (int i = 0; i <= p; ++i){
                for (int j = 0; j <= p; ++j){
                    tripletList.push_back(T(global_indices[i], global_indices[j], Element(i, j)));
                } 
            }
            

            std::vector<double> LinearElement=gauss_legendre_cuadrature_integral_linealForm(p,lowerLimit, upperLimit, nEvals, BasisFunsEvals, JacobianEvals, funcEvals, weights);
            
            if(PrintIntermediateValues){
                std::cout<< "--------------LinearElement----------------" << std::endl;
                for(unsigned int i=0; i<LinearElement.size(); i++){
                    std::cout << LinearElement[i] << ", ";
                }
                std::cout << std::endl;
            }
            for(unsigned int i=0; i<=p; i++){
                LinearForm(global_indices[i])+=LinearElement[i];
            }
            
        }
        
        global.setFromTriplets(tripletList.begin(), tripletList.end());

        if(PrintMatixes){
             std::cout<< "--------------GlobalMat----------------" << std::endl;
            for (int k = 0; k < global.outerSize(); ++k){
                for (Eigen::SparseMatrix<double>::InnerIterator it(global, k); it; ++it){
                    std::cout << "(" << it.row() << "," << it.col() << "): " << it.value() << "\n";
                }
            }
            std::cout<<Eigen::MatrixXd(global) <<std::endl;
            std::cout<< "--------------GlobalLinearForm----------------" << std::endl;
            for (int k = 0; k < LinearForm.size(); ++k){
                std::cout << LinearForm[k] << ", ";
            }
            std::cout << std::endl;
        }
       
        //Compute of the Inverse jacobians for boundary conditions on dirhclet/Robin conditions
        iMatrix<double> Aders, wders;
        iMatrix<double> Allders0 = CurveDerivsAlg1(n, p, KnotVector, NewCtrlPts, 0, 1);
        int auxInt=Allders0.GetNumCols()-1;
        Aders = Allders0.GetSubMat(0,1,0, auxInt-1);
        wders = Allders0.GetSubMat(0,1,auxInt, auxInt);
        iMatrix<double> AuxMat2 = RatCurveDerivs(Aders, wders, 1);
        std::vector<double> AuxVec= AuxMat2.GetRow(1);
        double JacobianEvals0=0;
        for(unsigned int i=0; i< AuxVec.size(); i++){
            JacobianEvals0+= AuxVec[i]*AuxVec[i];
        }

        
        double InverseJacobianEvals0=1/JacobianEvals0;


        iMatrix<double> Allders1 = CurveDerivsAlg1(n, p, KnotVector, NewCtrlPts, 1, 1);
        Aders = Allders1.GetSubMat(0,1,0, auxInt-1);
        wders = Allders1.GetSubMat(0,1,auxInt, auxInt);
        AuxMat2 = RatCurveDerivs(Aders, wders, 1);
        AuxVec= AuxMat2.GetRow(1);
        double JacobianEvals1=0;
        for(unsigned int i=0; i< AuxVec.size(); i++){
            JacobianEvals1+= AuxVec[i]*AuxVec[i];
        }
        double InverseJacobianEvals1=1/JacobianEvals1;




        //apply boundary conditions
        impose_Dirichlet_Condition(p, size-1, global, LinearForm, true, true);
        //impose_Newmann_Condition(size-1, LinearForm, true, true);
        double alpha = 3;
        double beta= 5;
        //impose_Robin_Condition(p, size-1, global, LinearForm, alpha, beta, true, false);

        if(PrintMatixes){
            std::cout<< "--------------GlobalMatPostCond----------------" << std::endl;
            for (int k = 0; k < global.outerSize(); ++k){
                for (Eigen::SparseMatrix<double>::InnerIterator it(global, k); it; ++it){
                    std::cout << "(" << it.row() << "," << it.col() << "): " << it.value() << "\n";
                }
            }
            std::cout<<Eigen::MatrixXd(global) <<std::endl;
        }
        
        if(PrintMatixes){
            std::cout<< "--------------GlobalLinearFormPostCond----------------" << std::endl;
            for (int k = 0; k < LinearForm.size(); ++k){
                std::cout << LinearForm[k] << ", ";
            }
            std::cout << std::endl;
        }
        

        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(global);
        if(solver.info() != Eigen::Success) {
            std::cout << "Error 1" << std::endl;
        }

        Eigen::VectorXd u = solver.solve(LinearForm);
        if(solver.info() != Eigen::Success) {
            std::cout << "Error 2" << std::endl;
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        std::cout << "Tiempo tomado: " << duration.count() << " nanosegundos\n";

        std::cout<< "--------------Solution----------------" << std::endl;
        for (int k = 0; k < u.size(); ++k){
            std::cout << u[k] << ", ";
        }
        std::cout << std::endl;


        std::string Text= "C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/SolEvals/EjSinNonConst_h=" + to_string(static_cast<int>(nElements)) +"_p="+ to_string(p)+".txt";
        writeFunction(u, KnotVector, NewCtrlPts.GetRow(2), n, p, lower_limit, upper_limit, Text);
        Text= "C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/SolEvals/EjSinNonConstDers_h=" + to_string(static_cast<int>(nElements)) +"_p="+ to_string(p)+".txt";
        writeFunctionDers(u, KnotVector,NewCtrlPts, NewCtrlPts.GetRow(2), n, p, lower_limit, upper_limit, Text);
        Text= "C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/SolEvals/EjSinNonConst_Analytic_h=" + to_string(static_cast<int>(nElements)) +".txt";
        writeFunction(analyticSol, lower_limit, upper_limit, nElements, NewCtrlPts, KnotVector, n, p, Text);
        Text= "C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/SolEvals/EjSinNonConstDers_Analytic_h=" + to_string(static_cast<int>(nElements)) +".txt";
        writeFunctionDers(analyticSolDers, lower_limit, upper_limit, nElements, NewCtrlPts, KnotVector, n, p, Text);
        */
        std::cout<<"fin 25"<<std::endl;
    }

    if(stoi(argv[1])==26){//test funciones auxiliares

        iMatrix<double> L= ReadDataFile_Matrix("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/EjPruebaFEM.txt");
        L.PrintMatrix();

        std::cout << L(0,1) << std::endl;
        
        double nElements = 100;

        std::vector<double> vec = createSubIntervals(nElements, 0,1);

        for(unsigned int i=0; i<=nElements; i++){
            std::cout << vec[i] << ", ";
        }
        std::cout<< std::endl;
        std::cout<<"fin 26"<<std::endl;
    }

    if(stoi(argv[1])==27){//test jacobian/ inverse jacobian

        int p=4;
        std::cout<< "----------------CtrlPts--------------" << std::endl;
        iMatrix<double> CtrlPts = ReadDataFile_Matrix("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/EjPruebaFEM.txt");        
        CtrlPts.PrintMatrix();
        std::vector<double> KnotVector={0,0,0,0,0,1,1,1,1,1};
        double lower_limit=KnotVector[0];
        double upper_limit=KnotVector[KnotVector.size()-1];
        std::vector<double> Weights={1,1,1,1,1};
        std::cout<< "----------WeightedCtrlPts--------------" << std::endl;
        iMatrix WeightedCtrlPts=WeightCtrlPts(CtrlPts,Weights);
        WeightedCtrlPts.PrintMatrix();
        
        int n= KnotVector.size()-p-2;
        std::cout<< "------------------------------" << std::endl;
        double EvalPoint = 0.4;
        std::cout << "EvalPoint = " << EvalPoint << std::endl;

        iMatrix<double> Aders, wders;
        iMatrix<double> Allders = CurveDerivsAlg1(n, p, KnotVector, WeightedCtrlPts, EvalPoint, 2);
        std::cout<< "----------Allders--------------" << std::endl;
        Allders.PrintMatrix();
        //std::cout << "i2" << std::endl;
        int auxInt=Allders.GetNumCols()-1;
        std::cout<< "----------Aders--------------" << std::endl;
        Aders = Allders.GetSubMat(0,1,0, auxInt-1);
        Aders.PrintMatrix();
        std::cout<< "----------wders--------------" << std::endl;
        wders = Allders.GetSubMat(0,1,auxInt, auxInt);
        wders.PrintMatrix();
        iMatrix<double> AuxMat2 = RatCurveDerivs(Aders, wders, 1);
        std::cout<< "----------RatCurveDerivs--------------" << std::endl;
        AuxMat2.PrintMatrix();
        std::vector<double> AuxVec= AuxMat2.GetRow(1);
        double AuxVal=0;
        for(unsigned int i=0; i< AuxVec.size(); i++){
            AuxVal+= AuxVec[i]*AuxVec[i];
        }
        AuxVal=sqrt(AuxVal);
        //std::cout << "i4" << std::endl;
        double JacobianEvals=AuxVal;
        double InverseJacobianEvals=1/AuxVal;

        std::cout<< "----------CurvePoint--------------" << std::endl;
        std::vector<double> CurvePoint=CurvePointRational(n,p,KnotVector, WeightedCtrlPts, EvalPoint);
        for(unsigned int i=0; i<CurvePoint.size(); i++){
            std::cout << CurvePoint[i] << ", ";
        }
        std::cout << std::endl;

        std::cout << "JacobianEval = " <<  std::setprecision(7) << JacobianEvals << std::endl;
        std::cout << "InverseJacobianEval = " << InverseJacobianEvals << std::endl;
        std::cout<<"fin 27"<<std::endl;
    }

    if(stoi(argv[1])==28){//test 1D problems in space

        bool PrintMatixes=true;
        bool PrintIntermediateValues=true;

        
        int p=4;
        double nElements = 16;
        //variables que no cambian en cada ejecución del bucle o que se declaran fuera para ir actualizandolas
        std::cout<< "----------------CtrlPts--------------" << std::endl;
        iMatrix<double> CtrlPts = ReadDataFile_Matrix("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/EjPruebaFEM.txt");        
        CtrlPts.PrintMatrix();
        std::vector<double> KnotVector={0,0,0,0,0,2,2,2,2,2};
        double lower_limit=KnotVector[0];
        double upper_limit=KnotVector[KnotVector.size()-1];
        std::vector<double> Weights={1,1,1,1,1};
        std::cout<< "----------WeightedCtrlPts--------------" << std::endl;
        iMatrix<double> WeightedCtrlPts=WeightCtrlPts(CtrlPts,Weights);
        WeightedCtrlPts.PrintMatrix();
        
        int n= KnotVector.size()-p-2;
        std::cout<< "---------------KnotVectorOld-----------------" << std::endl;
        for(unsigned int i=0; i<KnotVector.size(); i++){
            std::cout << KnotVector[i] << ", ";
        }
        std::cout << std::endl;
        std::cout << "n value:" << n << std::endl;
        
        std::vector<double> Insertions=createSubIntervals(nElements, lower_limit, upper_limit);
        Insertions=removeMatches(Insertions, KnotVector);
        iMatrix<double> NewCtrlPts=RefineKnotVectCurve(n,p,KnotVector, WeightedCtrlPts, Insertions, Insertions.size()-1);
        std::cout<< "----------------NewCtrlPtsW--------------" << std::endl;   
        NewCtrlPts.PrintMatrix();
        int size = KnotVector.size()-p-1;
        auto f = [](double x) { 
            double pi=3.141592653589793;
            return 10*10*pi*pi*sin(10*pi*x); 
        };

        auto analyticSol = [](double x){
            double pi=3.141592653589793;
            return sin(10*pi*x);
        };

        auto analyticSolDers = [](double x){
            double pi=3.141592653589793;
            return 10*pi*cos(10*pi*x);
        };

        n = KnotVector.size()-p-2;



        std::cout<< "---------------KnotVectorNew-----------------" << std::endl;
        for(unsigned int i=0; i<KnotVector.size(); i++){
            std::cout << KnotVector[i] << ", ";
        }
        double Lenght=IntegrateNormDer(lower_limit,upper_limit, 20, NewCtrlPts, KnotVector, n, p);

        std::cout << std::endl;
        std::cout << "n value:" << n << std::endl;
        
        Weights=NewCtrlPts.GetRow(NewCtrlPts.GetNumRows()-1);
        std::cout<< "---------------Weights-----------------" << std::endl;
        for(unsigned int i=0; i<Weights.size(); i++){
            std::cout << Weights[i] << ", ";
        }
        std::cout << std::endl;
        int dim=NewCtrlPts.GetNumRows()-1;
        int numElements= NewCtrlPts.GetNumCols();
        
        int nEvals=p+1;
        
        
        Eigen::SparseMatrix<double> global(size, size);
        Eigen::VectorXd LinearForm(size);
        LinearForm.setZero(); // Inicializa en cero
        
        
        std::cout<< "---------------Intervals-----------------" << std::endl;
        std::vector<double> Intervals=removeDuplicates(KnotVector);
        for(unsigned int i=0; i<Intervals.size(); i++){
            std::cout << Intervals[i] << ", ";
        }
        std::cout << std::endl;
        Eigen::VectorXd nodes, weights;
        legendre_pol(nEvals, nodes, weights);
            
        std::cout << "Nodos (raices):\n" << nodes.transpose() << "\n";
        std::cout << "Pesos:\n" << weights.transpose() << "\n"; 
        typedef Eigen::Triplet<double> T;
        std::vector<Eigen::Triplet<double>> tripletList;

        //iteraciones 
        auto start = std::chrono::high_resolution_clock::now();
        for(unsigned int index=0; index<Intervals.size()-1; index++){
            double lowerLimit=Intervals[index], upperLimit=Intervals[index+1];
            int span = FindSpan(KnotVector, lowerLimit,p,n);
            std::vector<int> global_indices(p+1);
            for(unsigned int w=0; w<=p; w++){
                global_indices[w]=span-p+w;
            }
            if(PrintIntermediateValues){

                std::cout << "global_indices: (";
                for(unsigned int w=0; w<p; w++){
                    std::cout <<global_indices[w] << ", ";
                }
                std::cout << global_indices[p] << ")" <<std::endl;
                std::cout << "Evaluation Interval: (" << Intervals[index] << ", " << Intervals[index+1]<< ")" <<std::endl;
                std::cout<< "--------------------------------------------" << std::endl;
            }
            
            iMatrix<double> BasisFunsEvals;
            iMatrix<double> DerBasisFunsEvals;
            std::vector<double> JacobianEvals;
            std::vector<double> InverseJacobianEvals;
            std::vector<double> funcEvals;
            
            D1_element_eval(n, span,p,nEvals, lowerLimit, upperLimit, KnotVector, Weights, nodes, NewCtrlPts, 
                            f, BasisFunsEvals, DerBasisFunsEvals, JacobianEvals, InverseJacobianEvals, funcEvals, Lenght);
            if(PrintIntermediateValues){
                std::cout<< "-------------BasisFunsEvals----------------" << std::endl;
                BasisFunsEvals.PrintMatrix();
                std::cout<< "-------------DerBasisFunsEvals----------------" << std::endl;
                DerBasisFunsEvals.PrintMatrix();
                std::cout<< "--------------JacobianEvals----------------" << std::endl;
                for(unsigned int i=0; i<JacobianEvals.size(); i++){
                    std::cout << JacobianEvals[i] << ", ";
                }
                std::cout << std::endl;
                std::cout<< "--------------InverseJacobianEvals----------------" << std::endl;
                for(unsigned int i=0; i<InverseJacobianEvals.size(); i++){
                    std::cout << InverseJacobianEvals[i] << ", ";
                }
                std::cout << std::endl;
                std::cout<< "--------------funcEvals----------------" << std::endl;
                for(unsigned int i=0; i<funcEvals.size(); i++){
                    std::cout << funcEvals[i] << ", ";
                }
                std::cout << std::endl;
                std::cout<< "--------------Element----------------" << std::endl;
            }

            iMatrix<double> Element=gauss_legendre_cuadrature_integral_bilinealForm(p,lowerLimit,upperLimit,nEvals, DerBasisFunsEvals, InverseJacobianEvals, weights);
            
            if(PrintIntermediateValues){
                Element.PrintMatrix();
                std::cout<< "-------------------------------------" << std::endl;
            }


            for (int i = 0; i <= p; ++i){
                for (int j = 0; j <= p; ++j){
                    tripletList.push_back(T(global_indices[i], global_indices[j], Element(i, j)));
                } 
            }
            

            std::vector<double> LinearElement=gauss_legendre_cuadrature_integral_linealForm(p,lowerLimit, upperLimit, nEvals, BasisFunsEvals, JacobianEvals, funcEvals, weights);
            
            if(PrintIntermediateValues){
                std::cout<< "--------------LinearElement----------------" << std::endl;
                for(unsigned int i=0; i<LinearElement.size(); i++){
                    std::cout << LinearElement[i] << ", ";
                }
                std::cout << std::endl;
            }
            for(unsigned int i=0; i<=p; i++){
                LinearForm(global_indices[i])+=LinearElement[i];
            }
            
        }
        
        global.setFromTriplets(tripletList.begin(), tripletList.end());

        if(PrintMatixes){
             std::cout<< "--------------GlobalMat----------------" << std::endl;
            for (int k = 0; k < global.outerSize(); ++k){
                for (Eigen::SparseMatrix<double>::InnerIterator it(global, k); it; ++it){
                    std::cout << "(" << it.row() << "," << it.col() << "): " << it.value() << "\n";
                }
            }
            std::cout<<Eigen::MatrixXd(global) <<std::endl;
            std::cout<< "--------------GlobalLinearForm----------------" << std::endl;
            for (int k = 0; k < LinearForm.size(); ++k){
                std::cout << LinearForm[k] << ", ";
            }
            std::cout << std::endl;
        }
       
        //Compute of the Inverse jacobians for boundary conditions on dirhclet/Robin conditions
        iMatrix<double> Aders, wders;
        iMatrix<double> Allders0 = CurveDerivsAlg1(n, p, KnotVector, NewCtrlPts, 0, 1);
        int auxInt=Allders0.GetNumCols()-1;
        Aders = Allders0.GetSubMat(0,1,0, auxInt-1);
        wders = Allders0.GetSubMat(0,1,auxInt, auxInt);
        iMatrix<double> AuxMat2 = RatCurveDerivs(Aders, wders, 1);
        std::vector<double> AuxVec= AuxMat2.GetRow(1);
        double JacobianEvals0=0;
        for(unsigned int i=0; i< AuxVec.size(); i++){
            JacobianEvals0+= AuxVec[i]*AuxVec[i];
        }

        
        double InverseJacobianEvals0=1/JacobianEvals0;


        iMatrix<double> Allders1 = CurveDerivsAlg1(n, p, KnotVector, NewCtrlPts, 1, 1);
        Aders = Allders1.GetSubMat(0,1,0, auxInt-1);
        wders = Allders1.GetSubMat(0,1,auxInt, auxInt);
        AuxMat2 = RatCurveDerivs(Aders, wders, 1);
        AuxVec= AuxMat2.GetRow(1);
        double JacobianEvals1=0;
        for(unsigned int i=0; i< AuxVec.size(); i++){
            JacobianEvals1+= AuxVec[i]*AuxVec[i];
        }
        double InverseJacobianEvals1=1/JacobianEvals1;




        //apply boundary conditions
        impose_Dirichlet_Condition(p, size-1, global, LinearForm, true, true);
        //impose_Newmann_Condition(size-1, LinearForm, true, true);
        double alpha = 3;
        double beta= 5;
        //impose_Robin_Condition(p, size-1, global, LinearForm, alpha, beta, true, false);

        if(PrintMatixes){
            std::cout<< "--------------GlobalMatPostCond----------------" << std::endl;
            for (int k = 0; k < global.outerSize(); ++k){
                for (Eigen::SparseMatrix<double>::InnerIterator it(global, k); it; ++it){
                    std::cout << setprecision(std::numeric_limits<double>::max_digits10) << "(" << it.row() << "," << it.col() << "): " << it.value() << "\n";
                }
            }
            std::cout<<Eigen::MatrixXd(global) <<std::endl;
        }
        
        if(PrintMatixes){
            std::cout<< "--------------GlobalLinearFormPostCond----------------" << std::endl;
            for (int k = 0; k < LinearForm.size(); ++k){
                std::cout  << setprecision(std::numeric_limits<double>::max_digits10) << LinearForm[k] << ", ";
            }
            std::cout << std::endl;
        }
        

        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(global);
        if(solver.info() != Eigen::Success) {
            std::cout << "Error 1" << std::endl;
        }

        Eigen::VectorXd u = solver.solve(LinearForm);
        if(solver.info() != Eigen::Success) {
            std::cout << "Error 2" << std::endl;
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        std::cout << "Tiempo tomado: " << duration.count() << " nanosegundos\n";

        std::cout<< "--------------Solution----------------" << std::endl;
        for (int k = 0; k < u.size(); ++k){
            std::cout << setprecision(std::numeric_limits<double>::max_digits10)<< u[k] << ", ";
        }
        std::cout << std::endl;

        //Pre-computations of necessary variables in order to write the functions
        nEvals=nElements;
        //nEvals =200;
        double auxeval1= (upper_limit - lower_limit)/2;
        double auxeval2= (upper_limit + lower_limit)/2;
        Eigen::VectorXd NewNodes, NewWeights;
        legendre_pol(nEvals, NewNodes, NewWeights);
        std::vector<double> time(nEvals);
        for(unsigned int i=0; i<nEvals; i++){
            time[i]=auxeval1*NewNodes(i)+auxeval2;
        }
        //std::cout<<"Longitud de la curva= "<< setprecision(std::numeric_limits<double>::max_digits10) <<  Lenght <<std::endl;
        std::string Text= "C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/SolEvals/EjSinNonConst"; //_h=" + to_string(static_cast<int>(nElements)) +"_p="+ to_string(p)+"_test2.txt"
        
        writeNumericFunction(u, KnotVector, Weights, n, p, nEvals, lower_limit, upper_limit, Text, NewNodes, NewWeights, time, nElements);
        //std::cout<<"Longitud de la curva= "<< setprecision(std::numeric_limits<double>::max_digits10) <<  Lenght <<std::endl;
        writeAnalyticFunction(analyticSol, analyticSolDers, lower_limit, upper_limit, nEvals, NewCtrlPts, KnotVector, n, p, Text, NewNodes, NewWeights, time, Lenght);
       
        /*
        writeFunction(analyticSol, lower_limit, upper_limit, 50, NewCtrlPts, KnotVector, n, p, Text, Lenght);
        Text= "C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/SolEvals/EjSinNonConstDers_Analytic_test2.txt";
        writeFunctionDers(analyticSolDers, lower_limit, upper_limit, 50, NewCtrlPts, KnotVector, n, p, Text, Lenght);
        */
        std::cout<<"Longitud de la curva= "<< setprecision(std::numeric_limits<double>::max_digits10) <<  Lenght <<std::endl;


        iMatrix<double> sol(dim,nEvals);
        std::vector<double> eval(dim);
        //std::cout<< "--------------------------------------------" << std::endl;
        //iMatrix WeightedCtrlPts=WeightCtrlPts(CtrlPts,Weights);
        //WeightedCtrlPts.PrintMatrix();
        for(unsigned int i = 0; i<nEvals; i++){
            //std::cout << i << std::endl;
            eval=CurvePointRational(KnotVector.size()-2-p,p,KnotVector,NewCtrlPts,time[i]);
            //std::cout << eval[0] << "|" << eval[1] << std::endl;
            sol.SetCol(i,eval);
        }
        Write_Curve(sol,"C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/Curvas/CurvaEjPruebaFEM.txt");

        std::cout<<"fin 28"<<std::endl;
    }

    if(stoi(argv[1])==29){//Fisical to parametric test
        int p=4;
        double nElements = 40;
        //variables que no cambian en cada ejecución del bucle o que se declaran fuera para ir actualizandolas
        std::cout<< "----------------CtrlPts--------------" << std::endl;
        iMatrix<double> CtrlPts = ReadDataFile_Matrix("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/EjPruebaFEM2.txt");        
        CtrlPts.PrintMatrix();
        std::vector<double> KnotVector={0,0,0,0,0,1,1,1,1,1};
        double lower_limit=KnotVector[0];
        double upper_limit=KnotVector[KnotVector.size()-1];
        std::vector<double> Weights={1,1,1,1,1};
        std::cout<< "----------WeightedCtrlPts--------------" << std::endl;
        iMatrix WeightedCtrlPts=WeightCtrlPts(CtrlPts,Weights);
        WeightedCtrlPts.PrintMatrix();
        
        int n= KnotVector.size()-p-2;
        std::cout<< "---------------KnotVectorOld-----------------" << std::endl;
        for(unsigned int i=0; i<KnotVector.size(); i++){
            std::cout << KnotVector[i] << ", ";
        }
        std::cout << std::endl;
        std::cout << "n value:" << n << std::endl;
        
        std::vector<double> Insertions=createSubIntervals(nElements,0,1);
        Insertions=removeMatches(Insertions, KnotVector);
        iMatrix<double> NewCtrlPts=RefineKnotVectCurve(n,p,KnotVector, WeightedCtrlPts, Insertions, Insertions.size()-1);
        std::cout<< "----------------NewCtrlPtsW--------------" << std::endl;   
        NewCtrlPts.PrintMatrix();
        int size = KnotVector.size()-p-1;
        n= KnotVector.size()-p-2;
        double evalPoint = 0.1;
        std::cout << "evalPoint= " << evalPoint << std::endl;
        double sol = Fisical_to_Parametric(evalPoint, lower_limit, upper_limit, 20, NewCtrlPts, KnotVector, n, p);
        std::cout << "sol= " << sol <<std::endl;
        std::vector<double> CurvePoint= CurvePointRational(n,p, KnotVector, NewCtrlPts, evalPoint);
        std::cout << "CurvePoint= ";
        for(unsigned int i=0; i<CurvePoint.size(); i++){
            std::cout << CurvePoint[i] << ", ";
        }
        std::cout <<std::endl;
    }
    if(stoi(argv[1])==30){//pre-D2 test
        
        //Initial definitions of variables
        int pY= 4;
        int pX= 4;
        int hY= 5;
        int hX= 7;
        std::vector<double> KnotVectorY = {0,0,0,0,0,1,1,1,1,1};
        std::vector<double> KnotVectorX = {0,0,0,0,0,0.5,1,1,1,1,1};

        auto analyticSol = [](double x, double y){
            double pi=3.141592653589793;
            return sin(pi*x)*sin(pi*y);
        };

        auto f = [](double x, double y){
            double pi=3.141592653589793;
            return 2*pi*pi*sin(pi*x)*sin(pi*y);
        };

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

        double lower_limitY=KnotVectorY[0];
        std::cout << "lower_limitY= " << lower_limitY << std::endl;
        double upper_limitY=KnotVectorY[KnotVectorY.size()-1];
        std::cout << "upper_limitY= " << upper_limitY << std::endl;
        double lower_limitX=KnotVectorX[0];
        std::cout << "lower_limit2= " << lower_limitX << std::endl;
        double upper_limitX=KnotVectorX[KnotVectorX.size()-1];
        std::cout << "upper_limit2= " << upper_limitX << std::endl;

        std::vector<double> data={1,1,1,1,1,1,
                                  1,1,1,1,1,1,
                                  1,1,1,1,1,1,
                                  1,1,1,1,1,1,
                                  1,1,1,1,1,1};

        iMatrix<double> Weights(5,6, data.data());
        std::vector<iMatrix<double>> CtrlPts = ReadDataFile_CtrlPts("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/EjSuperficieTriv.txt");
        
        std::cout<<"----CtrlPts----" << std::endl;
        for(unsigned int i=0; i<CtrlPts.size();i++){
            CtrlPts[i].PrintMatrix();
            std::cout<<"-------------------" << std::endl;
        }
        std::cout<<"---Matriz_Pesos---" << std::endl;
        Weights.PrintMatrix();
        std::cout<<"-----CtrlPtsW-----" << std::endl;
        std::vector<iMatrix<double>> CtrlPtsW=WeightCtrlPtsSurf(CtrlPts,Weights);
        for(unsigned int i=0; i<CtrlPtsW.size();i++){
            CtrlPtsW[i].PrintMatrix();
            std::cout<<"-------------------" << std::endl;
        }

        double N_Steps_Y = 100;
        double N_Steps_X = 100;

        double Step_Y=1/N_Steps_Y;
        double Step_X=1/N_Steps_X;

        std::vector<iMatrix<double>> eval(N_Steps_Y+1);
        for(unsigned int i=0; i<N_Steps_Y+1;i++){
            
            iMatrix<double> auxMatrix(3,N_Steps_X+1);
            for(unsigned int j=0; j< N_Steps_X+1;j++){
                
                std::vector<double> aux = SurfacePointRational(KnotVectorY.size()-2-pY,pY,KnotVectorY,KnotVectorX.size()-2-pX,pX,KnotVectorX,CtrlPtsW,i*Step_Y,j*Step_X);
                for(unsigned int k=0; k<3;k++){
                    auxMatrix(k,j)=aux[k];
            }
            }
            eval[i]=auxMatrix;
        }

        Write_BezierSurface(eval,"C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/Superficies/SuperficieRationalEjSuperficieTriv" );
        Write_CtrlPts(CtrlPts,"C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/EjSuperficieTriv");


        std::vector<double> InsertionsY=createSubIntervals(hY, lower_limitY, upper_limitY);
        InsertionsY=removeMatches(InsertionsY, KnotVectorY);
        std::cout<< "---------------InsertionsY-----------------" << std::endl;
        for(unsigned int i=0; i<InsertionsY.size(); i++){
            std::cout << InsertionsY[i] << ", ";
        }
        std::cout << std::endl;

        std::vector<double> InsertionsX=createSubIntervals(hX, lower_limitX, upper_limitX);
        InsertionsX=removeMatches(InsertionsX, KnotVectorX);
        std::cout<< "---------------InsertionsX-----------------" << std::endl;
        for(unsigned int i=0; i<InsertionsX.size(); i++){
            std::cout << InsertionsX[i] << ", ";
        }
        std::cout << std::endl;


        std::vector<iMatrix<double>> NewCtrlPtsWintermedieate= RefineKnotVectSurface(KnotVectorY.size()-2-pY,pY,KnotVectorY,KnotVectorX.size()-2-pX,pX,KnotVectorX, CtrlPtsW, true, InsertionsY, InsertionsY.size()-1);
        std::vector<iMatrix<double>> NewCtrlPtsW= RefineKnotVectSurface(KnotVectorY.size()-2-pY,pY,KnotVectorY,KnotVectorX.size()-2-pX,pX,KnotVectorX, NewCtrlPtsWintermedieate, false, InsertionsX, InsertionsX.size()-1);
        std::cout<<"------NewCtrlPtsW-----" << std::endl;
        for(unsigned int i=0; i<NewCtrlPtsW.size();i++){
            NewCtrlPtsW[i].PrintMatrix();
            std::cout<<"-------------------" << std::endl;
        }

        NewCtrlPtsWintermedieate.clear();
        std::cout<<"------NewWeights-----" << std::endl;
        iMatrix<double> NewWeights(NewCtrlPtsW.size(),NewCtrlPtsW[0].GetNumCols());
        for(unsigned int i=0; i<NewWeights.GetNumRows(); i++){
            for(unsigned int j=0; j<NewWeights.GetNumCols(); j++){
                NewWeights(i,j)=NewCtrlPtsW[i](3,j);
            }
        }
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


        std::cout<< "---------------Intervals1-----------------" << std::endl;
        std::vector<double> Intervals1=removeDuplicates(KnotVectorY);
        for(unsigned int i=0; i<Intervals1.size(); i++){
            std::cout << Intervals1[i] << ", ";
        }
        std::cout << std::endl;
        std::cout<< "---------------Intervals2-----------------" << std::endl;
        std::vector<double> Intervals2=removeDuplicates(KnotVectorX);
        for(unsigned int i=0; i<Intervals2.size(); i++){
            std::cout << Intervals2[i] << ", ";
        }
        std::cout << std::endl;

        std::vector<iMatrix<double>> evalRefinement(N_Steps_Y+1);
        for(unsigned int i=0; i<N_Steps_Y+1;i++){
            
            iMatrix<double> auxMatrix(3,N_Steps_X+1);
            for(unsigned int j=0; j< N_Steps_X+1;j++){
                
                std::vector<double> aux = SurfacePointRational(KnotVectorY.size()-2-pY,pY,KnotVectorY, KnotVectorX.size()-2-pX,pX,KnotVectorX,NewCtrlPtsW,i*Step_Y,j*Step_X);
                for(unsigned int k=0; k<3;k++){
                    auxMatrix(k,j)=aux[k];
            }
            }
            evalRefinement[i]=auxMatrix;
        }

        Write_BezierSurface(evalRefinement,"C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/Superficies/SuperficieRationalEjSuperficieTrivRefinement" );
        Write_CtrlPts(NewCtrlPtsW,"C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/EjSuperficieTrivRefinement");


        double tY= 0.25;
        double tX= 0.75;
        int MaxD=2;
        iMatrix<std::vector<double>> Allders = SurfaceDerivsAlg1(KnotVectorY.size()-2-pY,pY,KnotVectorY,KnotVectorX.size()-2-pX,pX,KnotVectorX,NewCtrlPtsW, tY, tX, MaxD);
        std::cout<<"--------Allders------" << std::endl;
        for(unsigned int i=0; i<=MaxD; i++){
            for(unsigned int j=0; j<=MaxD; j++){
                for(unsigned int k=0; k<Allders(0,0).size();k++){
                    std::cout << Allders(i,j)[k]<< " ";
                }
                std::cout<<"|";
            }
            std::cout << std::endl;
        }
        std::cout<<"--------Aders------" << std::endl;
        iMatrix<std::vector<double>> Aders(MaxD+1, MaxD+1);
        iMatrix<double> Wders(MaxD+1, MaxD+1);
        for(unsigned int i=0; i<=MaxD; i++){
            for(unsigned int j=0; j<=MaxD; j++){
                Aders(i,j)= std::vector<double>(Allders(i,j).begin(), Allders(i,j).begin()+3);
                Wders(i,j)=Allders(i,j)[3];
                for(unsigned int k=0; k<Aders(0,0).size();k++){
                    
                    std::cout << Aders(i,j)[k]<< " ";
                }
                std::cout<<"|";
            }
            std::cout << std::endl;
        }

        std::cout<<"--------Wders------" << std::endl;
        for(unsigned int i=0; i<=MaxD; i++){
            for(unsigned int j=0; j<=MaxD; j++){ 
                std::cout << Wders(i,j)<< " ";
                std::cout<<"|";
            }
            std::cout << std::endl;
        }


        /*
        std::cout<<"------Weights2-----" << std::endl;
        
        std::vector<iMatrix<double>> Weights2(Weights.GetNumRows());
        int dim=1;
        for(unsigned int i=0; i<CtrlPtsW.size();i++){
            iMatrix<double> aux(1, CtrlPtsW[0].GetNumCols());

                aux.SetRow(0,CtrlPtsW[i].GetRow(CtrlPtsW[0].GetNumRows()-1));
                
            Weights2[i]=aux;
            Weights2[i].PrintMatrix();
                    std::cout<<"-----------" << std::endl;
        }
        std::cout<<"------Wders-----" << std::endl;
        iMatrix<std::vector<double>> Wders = SurfaceDerivsAlg1(V1.size()-2-p1,p1,V1,V2.size()-2-p2,p2,V2,Weights2, t1, t2, MaxD);
        iMatrix<double> WdersAux(MaxD+1,MaxD+1);
        for(unsigned int i=0; i<=MaxD; i++){
            for(unsigned int j=0; j<=MaxD; j++){
                //for(unsigned int k=0; k<Wders(0,0).size();k++){
                    std::cout << Wders(i,j)[0]<< " ";
                    WdersAux(i,j)=Wders(i,j)[0];
               //}
                std::cout<<"|";
            }
            std::cout << std::endl;
        }
        WdersAux.PrintMatrix();
        */
        std::cout<< "---------------SurfRatDers-----------------" << std::endl;

        iMatrix<std::vector<double>> SurfRatDers= RatSurfaceDerivs(Aders, Wders, MaxD);
        std::cout<< "numRows: " <<SurfRatDers.GetNumRows()<< ", numCols: "<<SurfRatDers.GetNumCols()<< std::endl;
        for(unsigned int i=0; i<= MaxD; i++){
            for(unsigned int j=0; j<= MaxD; j++){
                for(unsigned int k=0; k<SurfRatDers(0,0).size();k++){
                    if(i+j<=2){
                         std::cout<<SurfRatDers(i,j)[k]<<" ";
                    }
                    else{
                        std::cout<<"0 ";
                    }
                   
                }
                std::cout<<"|";
            }
            std::cout<< std::endl;
        }
        
        int spanY=FindSpan(KnotVectorY, tY, pY, KnotVectorY.size()-pY-2);
        int spanX=FindSpan(KnotVectorX, tX, pX, KnotVectorX.size()-pX-2);
        std::cout << "spanY= " << spanY << std::endl;
        std::cout << "spanX= " << spanX << std::endl;
        

        std::vector<double> WeightsY=NewWeights.GetSubMat(0,NewWeights.GetNumRows()-1, spanY-pY, spanY-pY).GetCol(0);
        std::vector<double> WeightsX=NewWeights.GetSubMat(spanX-pX, spanX-pX,0,NewWeights.GetNumCols()-1).GetRow(0);

        std::cout<< "---------------WeightsY-----------------" << std::endl;
        for(unsigned int i=0; i<WeightsY.size(); i++){
            std::cout << WeightsY[i] << ", ";
        }
        std::cout << std::endl;

        std::cout<< "---------------WeightsX-----------------" << std::endl;
        for(unsigned int i=0; i<WeightsX.size(); i++){
            std::cout << WeightsX[i] << ", ";
        }
        std::cout << std::endl;

        std::vector<double> RatBasisY= RationalBasisFuns(spanY, tY, pY, KnotVectorY, WeightsY);
        std::vector<double> RatBasisX= RationalBasisFuns(spanX, tX, pX, KnotVectorX, WeightsX);

        std::cout<< "---------------RatBasisY-----------------" << std::endl;
        for(unsigned int i=0; i<RatBasisY.size(); i++){
            std::cout << RatBasisY[i] << ", ";
        }
        std::cout << std::endl;

        std::cout<< "---------------RatBasisX-----------------" << std::endl;
        for(unsigned int i=0; i<RatBasisX.size(); i++){
            std::cout << RatBasisX[i] << ", ";
        }
        std::cout << std::endl;

        iMatrix<double> DersBasisY = DerBasisFuns(spanY, tY, pY, 1, KnotVectorY);
        iMatrix<double> DersBasisX = DerBasisFuns(spanX, tX, pX, 1, KnotVectorX);
        iMatrix<double> activeWeights=NewWeights.GetSubMat(spanX-pX, spanX, spanY-pY, spanY);

        std::cout<< "---------------DersBasisY-----------------" << std::endl;
        DersBasisY.PrintMatrix();
        std::cout<< "---------------DersBasisX-----------------" << std::endl;
        DersBasisX.PrintMatrix();

        iMatrix<double> RationalizedDersBasisY=RationalizeDersBasisFuns(spanY, tY, pY, 1, KnotVectorY, DersBasisY, WeightsY);
        iMatrix<double> RationalizedDersBasisX=RationalizeDersBasisFuns(spanX, tX, pX, 1, KnotVectorX, DersBasisX, WeightsX);

        std::cout<< "---------------RationalizedDersBasisY-----------------" << std::endl;
        RationalizedDersBasisY.PrintMatrix();
        std::cout<< "---------------RationalizedDersBasisX-----------------" << std::endl;
        RationalizedDersBasisX.PrintMatrix();
        
        iMatrix<double> RatDersBasisY = RatDersBasisFuns(spanY, tY, pY, 1, KnotVectorY, WeightsY);
        iMatrix<double> RatDersBasisX = RatDersBasisFuns(spanX, tX, pX, 1, KnotVectorX, WeightsX);
        std::cout<< "---------------RatDersBasisY-----------------" << std::endl;
        RatDersBasisY.PrintMatrix();
        std::cout<< "---------------RatDersBasisX-----------------" << std::endl;
        RatDersBasisX.PrintMatrix();
        
        std::cout << "fin 30" << std::endl;
    }

    if(stoi(argv[1])==31){//D2 test
        
        bool PrintMatrixes=true;
        bool ShowValues=false;
        bool ShowIterations=false;
        bool ShowIterations2=false;

        //Initial definitions of variables
        int pY= 4;
        int pX= 4;
        int hY= 5;
        int hX= 5;
        int nEvalsY=pY+1;
        int nEvalsX=pX+1;
        int subMatSize=nEvalsX*nEvalsY;
        int nElementsY=hY+pY;
        int nElementsX=hX+pX;
        int size = nElementsY*nElementsX;
        
        Eigen::SparseMatrix<double> global(size, size);
        Eigen::VectorXd LinearForm(size);
        std::vector<double> KnotVectorY = {0,0,0,0,0,1,1,1,1,1};
        std::vector<double> KnotVectorX = {0,0,0,0,0,1,1,1,1,1};
                
        double lower_limitY=KnotVectorY[0];
        double upper_limitY=KnotVectorY[KnotVectorY.size()-1];
        double lower_limitX=KnotVectorX[0];
        double upper_limitX=KnotVectorX[KnotVectorX.size()-1];



        auto analyticSol = [](double x, double y){
            double pi=3.141592653589793;
            return sin(pi*x)*sin(pi*y);
        };

        auto f = [](double x, double y){
            double pi=3.141592653589793;
            return 2*pi*pi*sin(pi*x)*sin(pi*y);
        };

        std::vector<double> data={1,1,1,1,1,
                                  1,1,1,1,1,
                                  1,1,1,1,1,
                                  1,1,1,1,1,
                                  1,1,1,1,1};

        iMatrix<double> Weights(5,5, data.data());
        std::vector<iMatrix<double>> CtrlPts = ReadDataFile_CtrlPts("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/EjSuperficieTriv.txt");
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
        


        std::vector<iMatrix<double>> NewCtrlPtsWintermedieate= RefineKnotVectSurface(KnotVectorY.size()-2-pY,pY,KnotVectorY,KnotVectorX.size()-2-pX,pX,KnotVectorX, CtrlPtsW, true, InsertionsY, InsertionsY.size()-1);
        std::vector<iMatrix<double>> NewCtrlPtsW= RefineKnotVectSurface(KnotVectorY.size()-2-pY,pY,KnotVectorY,KnotVectorX.size()-2-pX,pX,KnotVectorX, NewCtrlPtsWintermedieate, false, InsertionsX, InsertionsX.size()-1);
        

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

        // Heurística de reserva por hilo para evitar muchas realocaciones.
        // Ajusta según memoria / experiencia.
        for (int t = 0; t < nthreads; ++t) {
            threadTriplets[t].reserve((hX * hY * nLoc * nLoc) / nthreads / 4 + 1000);
        }

        
        auto start = std::chrono::high_resolution_clock::now();
        //Precomputation of BasisFuns and its derivates.
        std::vector<int> SpanIndexY(IntervalsY.size()-1);
        std::vector<int> SpanIndexX(IntervalsX.size()-1);
        
        std::vector<std::vector<iMatrix<double>>> Basis_and_DersY=PreBasis_and_Ders(pY, 1, KnotVectorY, IntervalsY, nodesY,SpanIndexY);
        std::vector<std::vector<iMatrix<double>>> Basis_and_DersX=PreBasis_and_Ders(pX, 1, KnotVectorX, IntervalsX, nodesX,SpanIndexX);

        if(ShowValues){
            std::cout<< "---------------Basis_and_DersY-----------------" << std::endl;
            for(unsigned int i=0; i< Basis_and_DersY.size(); i++){
                std::cout<< "-  -  -  -  -  -  -  -  -  -  -  -" << std::endl;
                for(unsigned int j=0; j< nEvalsY; j++){
                    Basis_and_DersY[i][j].PrintMatrix();
                    std::cout<< "-  -  -  -  -  -  -  -  -  -  -  -" << std::endl;
                }
                std::cout<< "-+--+--+--+--+--+--+--+--+--+--+-" << std::endl;
            }
            std::cout<< "---------------Basis_and_DersX-----------------" << std::endl;
            for(unsigned int i=0; i< Basis_and_DersX.size(); i++){
                std::cout<< "-  -  -  -  -  -  -  -  -  -  -  -" << std::endl;
                for(unsigned int j=0; j< nEvalsX; j++){
                    Basis_and_DersX[i][j].PrintMatrix();
                    std::cout<< "-  -  -  -  -  -  -  -  -  -  -  -" << std::endl;
                }
                std::cout<< "-+--+--+--+--+--+--+--+--+--+--+-" << std::endl;
            }

        }
        //std::unordered_map<int, std::vector<std::pair<double,double>>> spanMapY=BuildSpanMap(IntervalsY, SpanIndexY);
        //std::unordered_map<int, std::vector<std::pair<double,double>>> spanMapX=BuildSpanMap(IntervalsX, SpanIndexX);

        // Double index matrix for easier and faster assembly
        
        #pragma omp parallel for collapse(2) schedule(dynamic)
        // Assembly of the matrix
        for (unsigned int i = 0; i < hX; i++) {
            for (unsigned int j = 0; j < hY; j++) {

                int tid = omp_get_thread_num();

                // Local buffers referidos por tid
                auto &localTriplets = threadTriplets[tid];
                auto &localLinear  = threadLinearForms[tid];

                iMatrix<double> K_e(subMatSize, subMatSize);
                if (ShowIterations) {
                    std::cout<<"--------Active Intervals---------"<< std::endl;
                    std::cout<< "IntervalX: (" << IntervalsX[i] << " " << IntervalsX[i+1]<<")"<< std::endl;
                    std::cout<< "IntervalY: (" << IntervalsY[j] << " " << IntervalsY[j+1]<<")"<< std::endl;
                    std::cout<< "SpanIndexX: (" << SpanIndexX[i] <<")"<< std::endl;
                    std::cout<< "SpanIndexY: (" << SpanIndexY[j] <<")"<< std::endl;
                }
                double auxeval1X= (IntervalsX[i+1] - IntervalsX[i])/2;
                double auxeval2X= (IntervalsX[i+1] + IntervalsX[i])/2;
                double auxeval1Y= (IntervalsY[i+1] - IntervalsY[i])/2;
                double auxeval2Y= (IntervalsY[i+1] + IntervalsY[i])/2;

                
                int spanU = SpanIndexX[i];
                int spanV = SpanIndexY[j];
                std::vector<double> WeightsY=NewWeights.GetSubMat(0,NewWeights.GetNumRows()-1, spanV-pY, spanV-pY).GetCol(0);
                std::vector<double> WeightsX=NewWeights.GetSubMat(spanU-pX, spanU-pX,0,NewWeights.GetNumCols()-1).GetRow(0);
                std::vector<iMatrix<double>> RationalizedDersBasisY=RationalizeDersBasisFunsVec(spanV, nEvalsY, nodesY, IntervalsY[j],IntervalsY[j+1], pY, 1, nY, KnotVectorY, Basis_and_DersY[j], WeightsY);
                std::vector<iMatrix<double>> RationalizedDersBasisX=RationalizeDersBasisFunsVec(spanU, nEvalsX,nodesX, IntervalsX[i], IntervalsX[i+1],pX, 1, nX ,KnotVectorX, Basis_and_DersX[i], WeightsX);
                if (ShowIterations2) {
                    std::cout<<"--------RationalizedDersBasisY---------"<< std::endl;
                    for(unsigned int iter=0; iter<nEvalsY; iter++){
                        RationalizedDersBasisY[iter].PrintMatrix();
                        std::cout<<"------------------------------"<< std::endl;
                    }
                    std::cout<<"--------RationalizedDersBasisX---------"<< std::endl;
                    for(unsigned int iter=0; iter<nEvalsX; iter++){
                        RationalizedDersBasisX[iter].PrintMatrix();
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
                        iMatrix<std::vector<double>> Allders = SurfaceDerivsAlg1(KnotVectorY.size()-2-pY,pY,KnotVectorY,KnotVectorX.size()-2-pX,pX,KnotVectorX,NewCtrlPtsW, EvalPointY, EvalPointX, 1);
                        
                        
                        iMatrix<std::vector<double>> Aders(2, 2);
                        iMatrix<double> Wders(2, 2);
                        for(unsigned int iteri=0; iteri<=1; iteri++){
                            for(unsigned int iterj=0; iterj<=1; iterj++){
                                Aders(iteri,iterj)= std::vector<double>(Allders(iteri,iterj).begin(), Allders(iteri,iterj).begin()+3);
                                Wders(iteri,iterj)=Allders(iteri,iterj)[3];
                            }
                        }
                        
                        //std::vector<double> SurfBasisRat=(RationalizedDersBasisX[a].GetRow(0)*RationalizedDersBasisY[b].GetRow(0)).Vectorize();
                        //std::vector<double> SurfBasisRatX=(RationalizedDersBasisX[a].GetRow(1)*RationalizedDersBasisY[b].GetRow(0)).Vectorize();
                        //std::vector<double> SurfBasisRatY=(RationalizedDersBasisX[a].GetRow(0)*RationalizedDersBasisY[b].GetRow(1)).Vectorize();
                        //iMatrix<std::vector<double>> SurfDers= RatSurfaceDerivs(Aders, Wders, 1);
                        
                        //iMatrix<double> PartialXX=SurfBasisRatX*SurfBasisRatX;
                        //iMatrix<double> PartialXY=SurfBasisRatX*SurfBasisRatY;
                        //iMatrix<double> PartialYY=SurfBasisRatY*SurfBasisRatY;
                        //double g11=dot_product(SurfDers(1,0),SurfDers(1,0));
                        //double g12=dot_product(SurfDers(1,0),SurfDers(0,1));
                        //double g22=dot_product(SurfDers(0,1),SurfDers(0,1));
                        //double detg=g11*g22-g12*g12;
                        //double Jinv=1/std::sqrt(detg);
                        //double auxconstant = weightsX[a]*weightsY[b]*Jinv;
                        //double g11c=  g22*auxconstant;
                        //double g12c= -g12*2*auxconstant;
                        //double g22c=  g11*auxconstant;
                        //iMatrix<double> Kq= g11c*PartialXX+g12c*PartialXY+g22c*PartialYY;
                        //K_e = K_e + Kq;
                        
                    }
                }

                // Emulation of local load vector (all ones)
                std::vector<double> F_e((pX+1)*(pY+1), 1.0);
                // Rellenar triplets en matriz global
                for(int a = 0; a <= pX; ++a){
                    
                    for(int b = 0; b <= pY; ++b){
                        int i_local = a*(pY+1) + b; 
                        int I_global = (spanU - pX + a) + (spanV - pY + b)*nElementsX;
                        localLinear[I_global] += F_e[i_local];
                        //std::cout<< "---------------------------"<< std::endl;
                        //std::cout << "I_global= " << I_global << std::endl;

                        LinearForm[I_global] += F_e[i_local];

                        for(int c = 0; c <= pX; ++c){
                            for(int d = 0; d <= pY; ++d){
                                
                                int j_local = c*(pY+1) + d; 
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

        std::cout << "fin 31" << std::endl;
    }
    return 0;
}


