/******************************************************************************* 
 *AUTHOR: Carlos Palomera Oliva
 *DATE: MARCH 2025
 *******************************************************************************/

// C:\Users\carlo\OneDrive\Escritorio\Uni\TFG\Programas\C++
// compilar con avx2 -mavx2
// g++ testeo.cpp iMatrix.cpp funciones.cpp -o testeo -mavx2
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
    if(stoi(argv[1])==1){ //Test CurvePointRational
        iMatrix CtrlPts = ReadDataFile_Matrix("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/Ej1.txt");        //cout<<"abierto"<<endl;
        CtrlPts.PrintMatrix();
        std::vector<float> KnotVector={0,0,0,0.2,0.4,0.6,0.8,1,1,1};
        std::vector<float> Weights={1,1,0.5,1,1,1,1};
        int p=2;
        int dim=2;
        float N_Steps = 1000;
        iMatrix<float> sol(dim,N_Steps+1);
        std::vector<float> eval(dim);
        float Step=1/(N_Steps);
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
        std::cout<<"fin 1"<<endl;
    }


    else if(stoi(argv[1])==2){ //test RationalBasisFuns

        int p=2;
        float N_Steps = 100;
        std::vector<float>KnotVector= {0,0,0,1,2,3,3,3};
        std::vector<float>WeightVector= {1,5,1,1,1};
        iMatrix<float> sol(KnotVector.size()-p-1,N_Steps+1);
        //std::vector<float> Nurb(NurbsVector.size());
        float Step=1/(N_Steps);
        for(unsigned int i = 0; i<N_Steps+1; i++){
            float t=Step*i*3;
            //std::cout << "Iteración iniciada " << i<< std::endl;
            int span = FindSpan(KnotVector, t, p, KnotVector.size()-2-p);
            std::vector<float> eval= RationalBasisFuns(span, t , p, KnotVector, WeightVector);
            for(int j = span-p; j<=span;j++){
                if(j>=0) {
                    sol.SetElement(j,i,eval[j-span+p]);
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
        Write_BasisFunctions(sol, "C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/RationalBasisFuncts");





        std::cout<<"fin 2"<<endl;


    }


    else if(stoi(argv[1])==3){ //test RationalBasisFuns

        int p=2;
        float N_Steps = 100;
        std::vector<float>KnotVector= {0,0,0,1,2,3,3,3};
        std::vector<float>WeightVector= {1,4,1,1,1};
        iMatrix<float> sol(KnotVector.size()-p-1,N_Steps+1);
        //std::vector<float> Nurb(NurbsVector.size());
        float Step=1/(N_Steps);
        for(unsigned int i = 0; i<N_Steps+1; i++){
            float t=Step*i*3;
            //std::cout << "Iteración iniciada " << i<< std::endl;
            int span = FindSpan(KnotVector, t, p, KnotVector.size()-2-p);
            std::vector<float> eval= BasisFuns(span, t , p, KnotVector);
            std::cout<< "j sera " << span-p<< std::endl;
            for(int j = span-p; j<=span;j++){
                if(j>=0) {
                    sol.SetElement(j,i,eval[j-span+p]);
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





        std::cout<<"fin 3"<<endl;


    }


     if(stoi(argv[1])==4){ //Test CurvePointRational
        iMatrix CtrlPts = ReadDataFile_Matrix("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/Ej1.txt");
        //cout<<"abierto"<<endl;
        CtrlPts.PrintMatrix();
                std::cout<< "--------------------------------------------" << std::endl;
        float t = 0.5;
        std::vector<float> KnotVector={0,0,0,0.2,0.4,0.6,0.8,1,1,1};
        std::vector<float> Weights={1,1,0.5,1,1,1,1};
        iMatrix<float> WeightsMat(1,Weights.size());
        WeightsMat.SetRow(0,Weights);
        WeightsMat.PrintMatrix();
    
        int p=2;
        int dim=2;
        float N_Steps = 100;
        iMatrix<float> sol(dim,N_Steps+1);
        std::vector<float> eval(dim);
        float Step=1/(N_Steps);
        std::cout<< "--------------------------------------------"<< std::endl;
        iMatrix WeightedCtrlPtsPreDrop=WeightCtrlPts(CtrlPts,Weights);
        iMatrix<float> WeightedCtrlPts(WeightedCtrlPtsPreDrop.GetNumRows()-1,WeightedCtrlPtsPreDrop.GetNumCols());
        for(unsigned int i=0; i< WeightedCtrlPts.GetNumRows();i++){
            WeightedCtrlPts.SetRow(i,WeightedCtrlPtsPreDrop.GetRow(i));
        }
        WeightedCtrlPts.PrintMatrix();
        iMatrix<float> Aders = CurveDerivsAlg1(KnotVector.size()-p-2, p, KnotVector, WeightedCtrlPts, t, 2);
        //std::cout << "p " << p<< " d " << 2 << std::endl;
        iMatrix<float> wders = CurveDerivsAlg1(KnotVector.size()-p-2, p, KnotVector, WeightsMat, t, 2);
        //std::cout << "p " << p<< " d " << 2 << std::endl;
        Aders.PrintMatrix();
        std::cout<<"----------------------------"<<std::endl;
        wders.PrintMatrix();
        std::cout<<"----------------------------"<<std::endl;


        iMatrix<float> CW= RatCurveDerivs(Aders, wders, 2);

        CW.PrintMatrix();


        std::cout<<"fin 4"<<endl;
    }

    if(stoi(argv[1])==5){ //Test SurfaceRational
        
        std::vector<iMatrix<float>> CtrlPts = ReadDataFile_CtrlPts("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/Ej2.txt");
        std::cout<<"----CtrlPts----" << std::endl;
        for(unsigned int i=0; i<CtrlPts.size();i++){
            CtrlPts[i].PrintMatrix();
            std::cout<<"-------------------" << std::endl;
        }
        iMatrix<float> Weights =ReadDataFile_Matrix("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Weights/Ej2Weights.txt");
        std::cout<<"---Matriz_Pesos---" << std::endl;
        Weights.PrintMatrix();
        std::cout<<"-----CtrlPtsW-----" << std::endl;
        std::vector<iMatrix<float>> CtrlPtsW=WeightCtrlPtsSurf(CtrlPts,Weights);
        for(unsigned int i=0; i<CtrlPtsW.size();i++){
            CtrlPtsW[i].PrintMatrix();
            std::cout<<"-------------------" << std::endl;
        }

        std::vector<float> V1={0,0,0,0.5,0.5,1,1,1};
        std::vector<float> V2={0,0,0,0,1,1,1,1};
        int p1=2;
        int p2=3;
        float N_Steps_1 = 100;
        float N_Steps_2 = 100;

        float Step_1=1/N_Steps_1;
        float Step_2=1/N_Steps_2;

        std::vector<iMatrix<float>> eval(N_Steps_1+1);
        for(unsigned int i=0; i<N_Steps_1+1;i++){
            
            iMatrix<float> auxMatrix(3,N_Steps_2+1);
            for(unsigned int j=0; j< N_Steps_2+1;j++){
                
                std::vector<float> aux = SurfacePointRational(V1.size()-2-p1,p1,V1, V2.size()-2-p2,p2,V2,CtrlPtsW,i*Step_1,j*Step_2);
                for(unsigned int k=0; k<3;k++){
                    auxMatrix.SetElement(k,j,aux[k]);
                }
            }
            eval[i]=auxMatrix;
        }

        Write_BezierSurface(eval,"C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/Superficies/SuperficieRationalEj2" );
        Write_CtrlPts(CtrlPts,"C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/Ej2");

    


        std::cout<<"fin 5"<<endl;
    }
    if(stoi(argv[1])==6){ //Test RatSurfaceDers
        
        float t1= 0.8f;
        float t2=0.8f;
        int p1=2;
        int p2=3;
        std::vector<float> V1={0,0,0,0.5,0.5,1,1,1};
        std::vector<float> V2={0,0,0,0,1,1,1,1};
        int MaxD=2;
        std::vector<iMatrix<float>> CtrlPts = ReadDataFile_CtrlPts("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/Ej2.txt");
        std::cout<<"----CtrlPts----" << std::endl;
        for(unsigned int i=0; i<CtrlPts.size();i++){
            CtrlPts[i].PrintMatrix();
            std::cout<<"-------------------" << std::endl;
        }
        iMatrix<float> Weights =ReadDataFile_Matrix("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Weights/Ej2Weights.txt");
        std::cout<<"---Matriz_Pesos---" << std::endl;
        Weights.PrintMatrix();
        std::cout<<"-----CtrlPtsW-----" << std::endl;
        std::vector<iMatrix<float>> CtrlPtsW=WeightCtrlPtsSurf(CtrlPts,Weights);
        
        for(unsigned int i=0; i<CtrlPtsW.size();i++){
            CtrlPtsW[i].PrintMatrix();
            std::cout<<"-------------------" << std::endl;
        }  
        
        std::vector<iMatrix<float>> CtrlPtsW2(CtrlPtsW.size());
        int dim=CtrlPtsW[0].GetNumRows()-1;
        for(unsigned int i=0; i<CtrlPtsW.size();i++){
            iMatrix<float> aux(dim, CtrlPtsW[0].GetNumCols());
            for(unsigned int j=0; j<dim;j++){
                aux.SetRow(j,CtrlPtsW[i].GetRow(j));
                
            }
            CtrlPtsW2[i]=aux;
        }



        iMatrix<std::vector<float>> Aders = SurfaceDerivsAlg1(V1.size()-2-p1,p1,V1,V2.size()-2-p2,p2,V2,CtrlPtsW2, t1, t2, MaxD);
        std::cout<<"--------Aders------" << std::endl;
        for(unsigned int i=0; i<=MaxD; i++){
            for(unsigned int j=0; j<=MaxD; j++){
                for(unsigned int k=0; k<Aders.GetElement(0,0).size();k++){
                    std::cout << Aders.GetElement(i,j)[k]<< " ";
                }
                std::cout<<"|";
            }
            std::cout << std::endl;
        }


        std::cout<<"------Weights2-----" << std::endl;
        
        std::vector<iMatrix<float>> Weights2(Weights.GetNumRows());
        dim=1;
        for(unsigned int i=0; i<CtrlPtsW.size();i++){
            iMatrix<float> aux(1, CtrlPtsW[0].GetNumCols());

                aux.SetRow(0,CtrlPtsW[i].GetRow(CtrlPtsW[0].GetNumRows()-1));
                
            Weights2[i]=aux;
            Weights2[i].PrintMatrix();
                    std::cout<<"-----------" << std::endl;
        }
        std::cout<<"------Wders-----" << std::endl;
        iMatrix<std::vector<float>> Wders = SurfaceDerivsAlg1(V1.size()-2-p1,p1,V1,V2.size()-2-p2,p2,V2,Weights2, t1, t2, MaxD);
        iMatrix<float> WdersAux(MaxD+1,MaxD+1);
        for(unsigned int i=0; i<=MaxD; i++){
            for(unsigned int j=0; j<=MaxD; j++){
                //for(unsigned int k=0; k<Wders.GetElement(0,0).size();k++){
                    std::cout << Wders.GetElement(i,j)[0]<< " ";
                    WdersAux.SetElement(i,j,Wders.GetElement(i,j)[0]);
                //}
                std::cout<<"|";
            }
            std::cout << std::endl;
        }
        WdersAux.PrintMatrix();

        iMatrix<std::vector<float>> SurfRatDers= RatSurfaceDerivs(Aders, WdersAux, MaxD);
        std::cout<< "numRows: " <<SurfRatDers.GetNumRows()<< ", numCols: "<<SurfRatDers.GetNumCols()<< std::endl;
        for(unsigned int i=0; i<= MaxD; i++){
            for(unsigned int j=0; j<= MaxD; j++){
                for(unsigned int k=0; k<SurfRatDers.GetElement(0,0).size();k++){
                    if(i+j<=2){
                         std::cout<<SurfRatDers.GetElement(i,j)[k]<<" ";
                    }
                    else{
                        std::cout<<"0 ";
                    }
                   
                }
                std::cout<<"|";
            }
            std::cout<< std::endl;
        }




        std::cout<<"fin 6"<<endl;    
    }

    if(stoi(argv[1])==7){ //Test Geometric Algorithms 1
        
        iMatrix CtrlPts = ReadDataFile_Matrix("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/Ej1.txt");
        //cout<<"abierto"<<endl;
        CtrlPts.PrintMatrix();
                std::cout<< "--------------------------------------------" << std::endl;
        std::vector<float> KnotVector={0,0,0,0.2,0.4,0.6,0.8,1,1,1};
        float element = 0.2f;
        int index=findInsertionIndex(KnotVector,element);
        int r = 1;
        std::vector<float> Weights={1,1,0.5,1,1,1,1};
        iMatrix<float> WeightsMat(1,Weights.size());
        WeightsMat.SetRow(0,Weights);
        WeightsMat.PrintMatrix();
        std::cout<< "--------------------------------------------"<< std::endl;
        iMatrix WeightedCtrlPts=WeightCtrlPts(CtrlPts,Weights);
        WeightedCtrlPts.PrintMatrix();
        std::cout<< "--------------------------------------------"<< std::endl;
        int p=2;
        int dim=2;
        iMatrix<float> NewCtrlPtsW= CurveKnotsIns(KnotVector.size()-2-p, p, KnotVector,WeightedCtrlPts, element, index, 0, r);
        NewCtrlPtsW.PrintMatrix();
        float N_Steps = 100;
        iMatrix<float> sol(dim,N_Steps+1);
        std::vector<float> eval(dim);
        float Step=1/(N_Steps);
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
        std::cout<< "fin 7";
    }
    if(stoi(argv[1])==8){ //Test Geometric Algorithms 1
        
        iMatrix CtrlPts = ReadDataFile_Matrix("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/Ej1.txt");
        std::vector<float> KnotVector={0,0,0,0.2,0.4,0.6,0.8,1,1,1};
        
        

        std::vector<float> Weights={1,1,0.5,1,1,1,1};
        iMatrix<float> WeightsMat(1,Weights.size());
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
        float N_Steps = 100;
        iMatrix<float> sol(dim,N_Steps+1);
        std::vector<float> eval(dim);
        float Step=1/(N_Steps);
        for(unsigned int i = 0; i<N_Steps+1; i++){
            //std::cout << i << std::endl;
            eval=CurvePntByCornerCut(KnotVector.size()-2-p,p,KnotVector,WeightedCtrlPts,i*Step);
            //std::cout << eval[0] << "|" << eval[1] << std::endl;
            sol.SetCol(i,eval);
        }
        Write_Curve(sol,"C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/Curvas/CurvaEj1CornerCut.txt");
        std::cout<<"fin 8";
    }

    if(stoi(argv[1])==9){
        std::vector<iMatrix<float>> CtrlPts = ReadDataFile_CtrlPts("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/Ej2.txt");
        std::cout<<"----CtrlPts----" << std::endl;
        for(unsigned int i=0; i<CtrlPts.size();i++){
            CtrlPts[i].PrintMatrix();
            std::cout<<"-------------------" << std::endl;
        }
        iMatrix<float> Weights =ReadDataFile_Matrix("C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Weights/Ej2Weights.txt");
        std::cout<<"---Matriz_Pesos---" << std::endl;
        Weights.PrintMatrix();
        std::cout<<"-----CtrlPtsW-----" << std::endl;
        std::vector<iMatrix<float>> CtrlPtsW=WeightCtrlPtsSurf(CtrlPts,Weights);
        for(unsigned int i=0; i<CtrlPtsW.size();i++){
            CtrlPtsW[i].PrintMatrix();
            std::cout<<"-------------------" << std::endl;
        }

        std::vector<float> V1={0,0,0,0.5,0.5,1,1,1};
        std::vector<float> V2={0,0,0,0,1,1,1,1};
        int p1=2;
        int p2=3;

        std::vector<iMatrix<float>> NewCtrlPtsW= SurfaceKnotIns(V1.size()-p1-2, p1, V1, V2.size()-p2-2,p2, V2, CtrlPtsW, true, 0.2f, findInsertionIndex(V2,0.2f) ,1,0);
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
        float N_Steps_1 = 100;
        float N_Steps_2 = 100;

        float Step_1=1/N_Steps_1;
        float Step_2=1/N_Steps_2;

        std::vector<iMatrix<float>> eval(N_Steps_1+1);
        for(unsigned int i=0; i<N_Steps_1+1;i++){
            
            iMatrix<float> auxMatrix(3,N_Steps_2+1);
            for(unsigned int j=0; j< N_Steps_2+1;j++){
                
                std::vector<float> aux = SurfacePointRational(V1.size()-2-p1,p1,V1, V2.size()-2-p2,p2,V2,NewCtrlPtsW,i*Step_1,j*Step_2);
                for(unsigned int k=0; k<3;k++){
                    auxMatrix.SetElement(k,j,aux[k]);
                }
            }
            eval[i]=auxMatrix;
        }
        std::cout<< std::endl;
        Write_BezierSurface(eval,"C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/Superficies/SuperficieRationalEj2Insertion" );
        std::vector<iMatrix<float>> NewCtrlPts = UnWeightCtrlPtsSurf(NewCtrlPtsW);
        std::cout<<"------NewCtrlPts-----" << std::endl;
        for(unsigned int i=0; i<NewCtrlPts.size();i++){
            NewCtrlPts[i].PrintMatrix();
            std::cout<<"-------------------" << std::endl;
        }
        Write_CtrlPts(NewCtrlPts,"C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/CtrlPts/Ej2Insertion");
        std::cout<<"fin 9";
    }
    return 0;
}