/******************************************************************************* 
 *AUTHOR: Carlos Palomera Oliva
 *DATE: MARCH 2025
 *******************************************************************************/

// C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\Programas\C++
// compilar con avx2 -mavx2
// g++ testeo.cpp iMatrix.cpp funciones.cpp -o testeo -mavx2
#include"SparseTensor.cpp" 
#include <iostream>
#include <vector>
#include"iMatrix.cpp"
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
            eval=CurveDerivsAlg1(Nurbs.size()-2-p,p,Nurbs,matrix, i*Step,dim);
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
            eval=CurveDerivsAlg2(Nurbs.size()-2-p,p,Nurbs,matrix, i*Step,dim);
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

        int p=2;
        double N_Steps = 100;
        std::vector<double> KnotVector={0,0,0,0.2,0.4,0.6,0.8,1,1,1};
        std::vector<double> WeightVector={1,1,0.5,1,1,1,1};
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
                    sol(j,i)=eval[j-span+p];                }
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
        iMatrix<double> Allders= CurveDerivsAlg1(KnotVector.size()-p-2, p, KnotVector, NewCtrlPts, t, 2);
        iMatrix<double> Aders = CurveDerivsAlg1(KnotVector.size()-p-2, p, KnotVector, WeightedCtrlPts, t, 2);
        //std::cout << "p " << p<< " d " << 2 << std::endl;
        iMatrix<double> wders = CurveDerivsAlg1(KnotVector.size()-p-2, p, KnotVector, WeightsMat, t, 2);
        //std::cout << "p " << p<< " d " << 2 << std::endl;
        std::cout<<"----------Allders----------"<<std::endl;
        Allders.PrintMatrix();
        std::cout<<"----------Aders----------"<<std::endl;
        Aders.PrintMatrix();
        std::cout<<"----------wders-------------"<<std::endl;
        wders.PrintMatrix();
        std::cout<<"-----------CW---------------"<<std::endl;


        iMatrix<double> CW= RatCurveDerivs(Aders, wders, 2);
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

    if(stoi(argv[1])==25){//test subMatrix1D

        bool PrintMatixes=false;
        bool PrintIntermediateValues=false;
        
        int p=4;
        double nElements = 80;
        //variables que no cambian en cada ejecución del bucle o que se declaran fuera para ir actualizandolas
        std::cout<< "----------------CtrlPts--------------" << std::endl;
        iMatrix<double> CtrlPts = ReadDataFile_Matrix("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/EjPruebaFEM.txt");        
        CtrlPts.PrintMatrix();
        std::vector<double> KnotVector={0,0,0,0,0,1,1,1,1,1};
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

        n = KnotVector.size()-p-2;

        std::cout<< "---------------KnotVectorNew-----------------" << std::endl;
        for(unsigned int i=0; i<KnotVector.size(); i++){
            std::cout << KnotVector[i] << ", ";
        }
        std::cout << std::endl;
        std::cout << "n value:" << n << std::endl;
        
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
        iMatrix<double> Allders0 = CurveDerivsAlg1(n, p, KnotVector, NewCtrlPts, 0, 2);
        int auxInt=Allders0.GetNumCols()-1;
        Aders = Allders0.GetSubMat(0,1,0, auxInt-1);
        wders = Allders0.GetSubMat(0,1,auxInt, auxInt);
        iMatrix<double> AuxMat2 = RatCurveDerivs(Aders, wders, 2);
        std::vector<double> AuxVec= AuxMat2.GetRow(1);
        double JacobianEvals0=0;
        for(unsigned int i=0; i< AuxVec.size(); i++){
            JacobianEvals0+= AuxVec[i]*AuxVec[i];
        }
        double InverseJacobianEvals0=1/JacobianEvals0;


        iMatrix<double> Allders1 = CurveDerivsAlg1(n, p, KnotVector, NewCtrlPts, 1, 2);
        Aders = Allders1.GetSubMat(0,1,0, auxInt-1);
        wders = Allders1.GetSubMat(0,1,auxInt, auxInt);
        iMatrix<double> AuxMat2 = RatCurveDerivs(Aders, wders, 2);
        std::vector<double> AuxVec= AuxMat2.GetRow(1);
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


        std::string Text= "C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/SolEvals/EjSin_h=" + to_string(static_cast<int>(nElements)) +"_p="+ to_string(p)+".txt";
        writeFunction(u, KnotVector, NewCtrlPts.GetRow(2), n, p, 0, 1, Text);
        Text= "C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/SolEvals/EjSin_Analytic_h=" + to_string(static_cast<int>(nElements)) +".txt";
        writeFunction(analyticSol, 0, 1, nElements, NewCtrlPts, KnotVector, n, p, Text);

        std::cout<<"fin 25"<<std::endl;
    }

    if(stoi(argv[1])==26){//test funciones auxiliares

        iMatrix<double> L= ReadDataFile_Matrix("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/EjPruebaFEM.txt");
        L.PrintMatrix();

        std::cout << L(0,1) << std::endl;
        
        double nElements = 100;

        std::vector<double> vec = createSubIntervals(nElements);

        for(unsigned int i=0; i<=nElements; i++){
            std::cout << vec[i] << ", ";
        }
        std::cout<< std::endl;
        std::cout<<"fin 26"<<std::endl;
    }
    return 0;
}


