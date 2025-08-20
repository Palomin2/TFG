/******************************************************************************* 
 *AUTHOR: Carlos Palomera Oliva
 *DATE: MARCH 2025
 *******************************************************************************/
#include "funciones.hpp"

/******************************************************************************* */

/******************************************************************************* */
iMatrix<double> ReadDataFile_Matrix(std::string direction){
    double Nrow,Ncol;
    std::string line;
    std::ifstream DataSet;

    DataSet.open(direction);
    if(DataSet.is_open()){

        if(std::getline(DataSet,line, ' ')){
            Nrow = stod(line);
        }
        if(std::getline(DataSet,line, '\n')){
            Ncol = stod(line);
        }

        iMatrix<double> matrix(Nrow,Ncol);


        for(unsigned int i=0; i<Nrow; i++){
            for(unsigned int j=0; j<Ncol-1;j++){ 
                if(std::getline(DataSet,line, ' ')){
                    matrix(i,j)=stod(line);
                }    
            }
            if(std::getline(DataSet,line, '\n')){
                matrix(i,Ncol-1)=stod(line);
            } 
        }
    DataSet.close();
	return matrix;
	}
    else{
        throw std::invalid_argument("Cannot open Data File.");
        iMatrix<double> m;
        return m;
    }

    
}
/******************************************************************************* */

/******************************************************************************* */
std::vector<iMatrix<double>> ReadDataFile_CtrlPts(std::string direction){
    double Ncol,NSplits;
    std::string line;
    std::ifstream DataSet;

    DataSet.open(direction);
    if(DataSet.is_open()){

        if(std::getline(DataSet,line, ' ')){
            Ncol = stod(line);
        }
        if(std::getline(DataSet,line, '\n')){
            NSplits = stod(line);
        }

        std::vector<iMatrix<double>> CtrlPts(NSplits);
        for(unsigned int k=0; k<NSplits;k++){

            CtrlPts[k].Resize(3,Ncol);

            for(unsigned int i=0; i<3; i++){

                for(unsigned int j=0; j<Ncol-1;j++){

                    if(std::getline(DataSet,line, ' ')){

                        CtrlPts[k](i,j)=stod(line);
                    }    
                }
                if(std::getline(DataSet,line, '\n')){
                    CtrlPts[k](i,Ncol-1)=stod(line);
                } 
            }
            getline(DataSet,line,'\n');

        }

        
    DataSet.close();
	return CtrlPts;
	}
    else{
        throw std::invalid_argument("Cannot open Data File.");
        std::vector<iMatrix<double>> T;
        return T;
    }
}


/******************************************************************************* */

/******************************************************************************* */
std::vector<double> BernsteinPolinomial(double t, int n){
    std::vector<double> Bnew(n);
    std::vector<double> Bold(n);
    std::fill(Bold.begin(), Bold.end(), 0);
    Bold[0]=1; 
    double t1=1-t; 
    for(unsigned int j=1;j<n;j++){

        std::fill(Bnew.begin(), Bnew.end(), 0);

        for(unsigned int i=0; i<j;i++){

            Bnew[i]=Bnew[i]+t1*Bold[i];
            Bnew[i+1]=Bnew[i+1]+t*Bold[i];
        }
        
        Bold=Bnew;

    }
    return Bnew;

}

/******************************************************************************* */

/******************************************************************************* */
std::vector<double> BezierCurve_Point_t(iMatrix<double> &PtosControl, double t){
    int dim = PtosControl.GetNumRows();
    std::vector<double> Curve(dim);        
    int n=PtosControl.GetNumCols();
    std::vector<double> B = BernsteinPolinomial(t,n);
    for(unsigned int j=0; j<dim; j++){
        for(unsigned int i=0; i<n;i++){
            Curve[j]=Curve[j]+B[i]*PtosControl(j,i);
        }
    }
    return Curve;
}

/******************************************************************************* */

/******************************************************************************* */
iMatrix<double> BezierCurve(iMatrix<double> &PtosControl, double N_Steps){
    
    iMatrix<double> Curve(PtosControl.GetNumRows(),N_Steps+1);
    double Step=1/(N_Steps);
    for(unsigned int i=0; i<N_Steps+1;i++){
        Curve.SetCol(i,BezierCurve_Point_t(PtosControl,i*Step));
    }
    return Curve;
}

/******************************************************************************* */

/******************************************************************************* */
std::vector<double> BezierSurface_Point_t1_t2(std::vector<iMatrix<double>> &PtosControl, double t1, double t2){
    int dim = 3;
    std::vector<double> Curve(dim);

    int N_Splits = PtosControl.size();
    int n=PtosControl[0].GetNumCols();
    std::vector<double> B1 = BernsteinPolinomial(t1,N_Splits);
    std::vector<double> B2 = BernsteinPolinomial(t2,n);

    for(unsigned int k=0; k<N_Splits;k++){
        for(unsigned int j=0; j<dim; j++){
            for(unsigned int i=0; i<n;i++){
                Curve[j]=Curve[j]+B1[k]*B2[i]*PtosControl[k](j,i);
            }
        }
    }
    
    return Curve;
}

/******************************************************************************* */

/******************************************************************************* */
std::vector<iMatrix<double>> BezierSurface(std::vector<iMatrix<double>> &PtosControl, double N_Steps){
    std::vector<iMatrix<double>> Surface(N_Steps+1);
    double Step=1/(N_Steps);

    for(unsigned int k=0; k<N_Steps+1;k++){
        Surface[k].Resize(3,N_Steps+1);
        for(unsigned int i=0; i<N_Steps+1;i++){
            Surface[k].SetCol(i,BezierSurface_Point_t1_t2(PtosControl,i*Step,k*Step));
        }
    }
    return Surface;
    
}

/******************************************************************************* */

/******************************************************************************* */
void Write_BezierSurface(std::vector<iMatrix<double>> &Tensor, std::string name){
    std::ofstream fichero;
    char coord[3]={'X','Y','Z'};
    int iterj = Tensor[0].GetNumCols();
    int iterk = Tensor.size();

    for(unsigned int i=0; i<3;i++){ 
        fichero.open(name+coord[i]+".txt");
        for(unsigned int j=0; j<iterj;j++){
            for(unsigned int k=0; k<iterk;k++){
                fichero<<Tensor[k](i,j)<<" ";
            }
            fichero<<std::endl;
        }
        fichero.close();
        
    }
}

/******************************************************************************* */

/******************************************************************************* */
void Write_CtrlPts(std::vector<iMatrix<double>> &CtrlPts, std::string name){
    std::ofstream fichero;
    char coord[3]={'X','Y','Z'};
    int iterj = CtrlPts[0].GetNumCols();
    int iterk = CtrlPts.size();

    for(unsigned int i=0; i<3;i++){ 
        fichero.open(name+coord[i]+".txt");
        for(unsigned int j=0; j<iterj;j++){
            for(unsigned int k=0; k<iterk;k++){
                fichero<<CtrlPts[k](i,j)<<" ";
            }
            fichero<<std::endl;
        }
        fichero.close();
        
    }

}



/************************************************************************************
*
* Classic funtions from the Nurbs Book 
*
************************************************************************************/

int FindSpan(std::vector<double> &KnotVector, double t, int p, int n){

    if(KnotVector[n+1]==t){
        return n;
    }

    else{
        int low =p;
        int high=n+1;
        int mid=(low+high)/2;
        while(t<KnotVector[mid] || t>=KnotVector[mid+1]){
            if(t<double(KnotVector[mid])){
                //std::cout<<"Condición t<U[mid]"<< "t="<<t<< " mid="<<mid<< ", U[mid]="<< KnotVector[mid]<<std::endl;
                high=mid;
            }
            else{
                //std::cout<<"Condición t>=U[mid]"<< "t="<<t<< " mid="<<mid<< ", U[mid]="<< KnotVector[mid]<<std::endl;
                low=mid;
            }
            mid=(low+high)/2;
        }
        return mid;
    }
}

/******************************************************************************* */

/******************************************************************************* */
std::vector<double> BasisFuns(int i, double t, int p, std::vector<double> &KnotVector){
    std::vector<double> sol(p+1);
    std::vector<double> left(p+1);
    std::vector<double> right(p+1);
    double temp,saved;
    //std::cout<<"Entramos en basisFuns"<<std::endl;
    sol[0]=1.0f;
    for(unsigned int j=1; j<=p; j++){
        //std::cout<< "j: " << j<< std::endl;
        left[j]=t-KnotVector[i+1-j];
        right[j]=KnotVector[i+j]-t;
        saved = 0.0f;
        for(unsigned int r=0; r<j; r++){
            //std::cout<< "r: " << r<< std::endl;
            temp= sol[r]/(right[r+1]+left[j-r]);
            sol[r]=saved+right[r+1]*temp;
            saved=left[j-r]*temp;
        }
        sol[j]=saved;
    }
    //std::cout<<"Salimos de basisFuns"<<std::endl;
    return sol;
}

/******************************************************************************* */

/******************************************************************************* */



iMatrix<double> AllBasisFuns(int i, double t, int p, std::vector<double> &KnotVector){


    iMatrix<double> sol(p+1,p+1);
    std::vector<double> left(p+1);
    std::vector<double> right(p+1);
    sol(0,0)=1.0f;
    double temp,saved;
  
    for(unsigned int j=1; j<=p; j++){
        
        left[j]=t-KnotVector[i+1-j];
        right[j]=KnotVector[i+j]-t;
        saved = 0.0f;
        for(unsigned int r=0; r<j; r++){
            
            temp= sol(j-1,r)/(right[r+1]+left[j-r]);
            sol(j,r)=saved+right[r+1]*temp;
            saved=left[j-r]*temp;
        }
        sol(j,j)=saved;
        
    }
    
    return sol;
}

/******************************************************************************* */

/******************************************************************************* */
iMatrix<double> DerBasisFuns(int span, double t, int p, int k,std::vector<double> &KnotVector){
    iMatrix<double> ders(p+1,p+1);
    iMatrix<double> ndu(p+1,p+1);
    double temp,saved;
    ndu(0,0)=1.0f;
    // Almacenamiento para diferencias izquierda y derecha
    std::vector<double> left(p + 1), right(p + 1);

    // Llenar la tabla de funciones base (algoritmo recursivo)
    for (int j = 1; j <= p; ++j) {
        left[j] = t - KnotVector[span + 1 - j];
        right[j] = KnotVector[span + j] - t;
        double saved = 0.0;

        for (int r = 0; r < j; ++r) {
            ndu(j,r)= right[r + 1] + left[j - r];
            double temp= ndu(r,j-1)/ndu(j,r);
            ndu(r,j)=saved + right[r + 1] * temp ;
            saved = left[j - r] * temp;
        }
        ndu(j,j)=saved;
    }

    // Inicializar derivadas (0-ésima derivada = función base)
    for (int j = 0; j <= p; ++j) {
        ders(0,j)=ndu(j,p);
    }

    // Calcular derivadas de orden superior
    for (int r = 0; r <= p; ++r) {
        int s1 = 0, s2 = 1; // Alternar filas en a
        iMatrix<double> a(2,p+1);
        a(0,0)=1.0f;

        for (int i = 1; i <= k; ++i) {
            double d = 0.0;
            int rk = r - i, pk = p - i;
            if (r >= i) {
                a(s2,0)= a(s1,0)/ndu(pk+1,rk);
                d= a(s2,0)*ndu(rk,pk);
            }
            int j1 = (rk >= -1) ? 1 : -rk;
            int j2 = (r - 1 <= pk) ? i - 1 : p - r;

            for (int j = j1; j <= j2; ++j) {
                a(s2,j)= (a(s1,j)-a(s1,j-1))/ndu(pk+1,rk+1)   ;
                d += a(s2,j)*ndu(rk+j,pk);
            }
            if (r <= pk) {
                a(s2,i)= -a(s1,i-1)/ndu(pk+1,r)     ;
                d += a(s2,i)*ndu(r,pk);
            }
            ders(i,r)=d;
            std::swap(s1, s2); // Alternar filas
        }
        
    }

    // Ajustar multiplicadores
    int r = p;
    for (int i = 1; i <= k; ++i) {
        for (int j = 0; j <= p; ++j) {
            ders(i,j)= ders(i,j)*r;
        }
        r *= (p - i);
    }

    return ders;
}

/******************************************************************************* */

/******************************************************************************* */
void Write_BasisFunctions(iMatrix<double> &BasisFunctions, std::string name){
    std::ofstream fichero;
    fichero.open(name+".txt");
    for(unsigned int j=0; j<BasisFunctions.GetNumRows();j++){
        for(unsigned int k=0; k<BasisFunctions.GetNumCols();k++){
            fichero<<BasisFunctions(j,k)<<" ";
        }
        fichero<<std::endl;
    }
    fichero.close();
}

/************************************************************************************
*
* Classic funtions for p spline curves/surfaces
*
************************************************************************************/

std::vector<double> CurvePoint(int n, int p, std::vector<double> &KnotVector, iMatrix<double> &CtrlPts, double t){
    int span =FindSpan(KnotVector, t, p, n);
    std::vector<double> Funs = BasisFuns(span,t,p,KnotVector);
    
    int itermax=CtrlPts.GetNumRows();
    std::vector<double> C(itermax);
    for(unsigned int i = 0; i<=p;i++){
        double aux=Funs[i];
        
        for(unsigned int j=0;j<itermax;j++){
            
            C[j]=C[j]+ aux*(CtrlPts(j,span-p+i));
            
        }  
       
    }

    return C;
}

/******************************************************************************* */

/******************************************************************************* */
iMatrix<double> CurveDerivsAlg1(int n, int p,std::vector<double> &KnotVector, iMatrix<double> &CtrlPts, double t, int d){
    int du= std::min(d,p);
    int span =FindSpan(KnotVector, t, p, n);
    int itermax=CtrlPts.GetNumRows();
    //std::cout << "du " << du<< std::endl;
    iMatrix<double> CK(du, itermax); //inicializa a 0 no hace falta dar valores
    iMatrix<double> nders= DerBasisFuns(span, t, p, du,KnotVector);
    for(unsigned int k=0; k<=du;k++){
        for(unsigned int i=0; i<=p;i++){
            double aux =nders(k,i);

            //CK[k] = CK[k] + nders[k] [j]*P[span-p+j]
            for(unsigned int j=0;j<itermax;j++){
            
                CK(k,j)=CK(k,j)+aux*(CtrlPts(j,span-p+i)) ;      
                
            } 
        }
    }
    return CK;

}

/******************************************************************************* */

/******************************************************************************* */
iMatrix<std::vector<double>> CurveDerivsCpts(int n, int p,std::vector<double> &KnotVector, iMatrix<double> &CtrlPts,  int d, int r1, int r2){
    int r= r2-r1;
    iMatrix<std::vector<double>> sol(d+1,r+1);
    for(unsigned int i=0;i<=r;i++){
        sol(0,i)= CtrlPts.GetCol(r1+i);
    }
    for(unsigned int k=1;k<=d;k++){
        int tmp = p-k+1;
        for(unsigned int i=0;i<=r-k;i++){
            
            int iter=CtrlPts.GetNumCols(); //dimension del espacio
            std::vector<double> aux(iter);
            double denominator=KnotVector[r1+i+p+1]-KnotVector[r1+i+k];

            for(unsigned int j=0; j<iter; j++){
                aux[j]=tmp*(sol(k-1,i+1)[j]-sol(k-1,i)[j])/denominator;
            }

            sol(k,i)= aux  ;
        }
    }

    return sol;
}

/******************************************************************************* */

/******************************************************************************* */
iMatrix<double> CurveDerivsAlg2(int n, int p,std::vector<double> &KnotVector, iMatrix<double> &CtrlPts, double t, int d){ //not working atm
    int du= std::min(d,p);
    int span =FindSpan(KnotVector, t, p, n);
    int itermax=CtrlPts.GetNumRows(); //dimension del espacio
    iMatrix<double> CK(itermax, du); //inicializa a 0 no hace falta dar valores
    iMatrix<double> nders= AllBasisFuns(span, t, p, KnotVector);
    iMatrix<std::vector<double>> DervCtps=CurveDerivsCpts(n,p,KnotVector,CtrlPts,du,span-p,span);
    for(unsigned int k=0; k<=du;k++){
        for(unsigned int i=0; i<=p-k;i++){
            double aux =nders(i,p-k);
            for(unsigned int j=0;j<itermax;j++){
            
                CK(k,j)=CK(k,j)+aux*DervCtps(k,i)[j] ;      
                
            } 
        }
    }
    return CK;

}

/******************************************************************************* */

/******************************************************************************* */
void Write_Curve(iMatrix<double> &Curve, std::string filename){
    std::ofstream fichero;
    int iterMaxi=Curve.GetNumCols();
    int iterMaxj=Curve.GetNumRows();
    fichero.open(filename);
    for(unsigned int j=0; j<iterMaxj;j++){
        for(unsigned int i=0; i<iterMaxi;i++){
            fichero<<Curve(j,i)<<" ";
        }
        fichero<< std::endl;
    }
    fichero.close();
}

/******************************************************************************* */

/******************************************************************************* */
std::vector<double> SurfacePoint(int n1, int p1, std::vector<double> &KnotVector1, int n2, int p2,  std::vector<double> &KnotVector2, std::vector<iMatrix<double>> &CtrlPts, double t1, double t2){
    int span1=FindSpan(KnotVector1,t1,p1,n1);
    int span2=FindSpan(KnotVector2,t2,p2,n2);
    std::vector<double> Basis1=BasisFuns(span1,t1,p1,KnotVector1);
    std::vector<double> Basis2=BasisFuns(span2,t2,p2,KnotVector2);
    int ind1 = span1-p1;
    std::vector<double> Pt= {0,0,0};
    for(unsigned int l=0; l<=p2; l++){

        std::vector<double> temp= {0,0,0};
        
        int ind2= span2 - p2+l;
        for(unsigned int k=0;k<=p1;k++){
            double auxBasis1= Basis1[k];
            //std::vector<double> auxCtrlPt = CtrlPts(ind1+k,ind2); 
            for(unsigned int i=0; i<3;i++){
                //std::cout<<"i="<<i<<",k="<<k<<",l="<<l<<std::endl;
                temp[i]=temp[i]+auxBasis1*CtrlPts[ind1+k](i,ind2);
            }
        }

        double auxBasis2=Basis2[l];
        for(unsigned int i=0; i<3;i++){
            Pt[i]=Pt[i]+auxBasis2*temp[i];
        }
    }
    return Pt;
}

/******************************************************************************* */

/******************************************************************************* */
iMatrix<std::vector<double>> SurfaceDerivsAlg1(int n1, int p1, std::vector<double> &KnotVector1, int n2, int p2,  std::vector<double> &KnotVector2, std::vector<iMatrix<double>> &CtrlPts, double t1, double t2, int d){
    unsigned int du=std::min(d,p1);
    unsigned int dv=std::min(d,p2);

    unsigned int dim = CtrlPts[0].GetNumRows();
    //std::cout<< "dim: "<< dim << std::endl;
    iMatrix<std::vector<double>> Sol(du+1,dv+1); 
    const std::vector<double> zero(dim);
    for(unsigned int k=0; k<=du; k++){
        for(unsigned int j=0; j<=dv;j++){
            Sol(k,j)= zero ;
        }
    }
    //std::cout<< "2" << std::endl;
    int span1=FindSpan(KnotVector1,t1,p1,n1);
    int span2=FindSpan(KnotVector2,t2,p2,n2);
    
    iMatrix<double> Ders1 = DerBasisFuns(span1, t1, p1, du, KnotVector1);
    iMatrix<double> Ders2 = DerBasisFuns(span2, t2, p2, dv, KnotVector2);
    iMatrix<double> temp(dim,p2+1); 
    //std::cout<< "3" << std::endl;
    for(unsigned int k=0; k<=du; k++){
            //std::cout<< "4" << std::endl;
        
        for(unsigned int s=0;s<=p2;s++){
                //std::cout << "s: "<< s << std::endl;
            temp.SetCol(s,zero);
            for(unsigned int r=0; r<=p1; r++){
                //    std::cout<< "6" << std::endl;
                //temp[s] = temp[s] + Nu[k] [r]*P[uspan-p+r] [vspan-q+s];
                double aux= Ders1(k,r);
                for(unsigned int j=0; j<dim;j++){
                    //    std::cout<< "7" << std::endl;
                    //std::cout<< CtrlPts[span1-p1+r](j,span2 -p2+s) << "|";
                    temp(j,s)= temp(j,s)+aux*CtrlPts[span1-p1+r](j,span2 -p2+s);
                }
                

            }
            //std::cout << "---------------------" << std::endl;
        }
        //std::cout << "8" << std::endl;
        int dd=std::min(d-k,dv);
        for(unsigned int l=0; l<= dd; l++){
            Sol(k,l)=zero;
            for(unsigned int s=0;s<=p2;s++){
                double aux = Ders2(l,s);
                std::vector<double> auxSol= Sol(k,l);
                std::vector<double> aux2(dim);
                for(unsigned int j=0; j<dim; j++){
                    aux2[j]=auxSol[j]+aux*temp(j,s);
                }
                Sol(k,l)= aux2;

            }
        }
    }

    return Sol;
}

/************************************************************************************
*
* Classic and support funtions for NURBS curves/surfaces
*
************************************************************************************/
iMatrix<double> WeightCtrlPts(const iMatrix<double> &CtrlPts,const std::vector<double> &Weights){
    int dim= CtrlPts.GetNumRows()+1;
    int size = CtrlPts.GetNumCols();
    iMatrix<double> CtrlPtsWeighted(dim,size);
    for(unsigned int i=0; i<size;i++){
        for(unsigned int j=0; j<dim-1; j++){
            CtrlPtsWeighted(j,i)=CtrlPts(j,i)*Weights[i];
        }
        CtrlPtsWeighted(dim-1,i)=Weights[i];
    }

    return CtrlPtsWeighted;
}

/******************************************************************************* */

/******************************************************************************* */

iMatrix<double> UnWeightCtrlPts(iMatrix<double> &CtrlPtsW){
    int NumPts = CtrlPtsW.GetNumCols();
    int dim = CtrlPtsW.GetNumRows()-1;
    iMatrix<double> CtrlPts(dim,NumPts);
    for(unsigned int j=0; j< NumPts;j++){
        for(unsigned int i=0; i<dim; i++){
            if(CtrlPtsW(dim,j)==0){
                //nothing, the declararion automatically sets it to 0
            }
            else{
                CtrlPts(i,j)=CtrlPtsW(i,j)/CtrlPtsW(dim,j);
            }
        }
    }
    return CtrlPts;
}

/******************************************************************************* */

/******************************************************************************* */
std::vector<iMatrix<double>> WeightCtrlPtsSurf(const std::vector<iMatrix<double>> &CtrlPts,const iMatrix<double> &Weights){
    int dim = CtrlPts.size();
    std::vector<iMatrix<double>> WeightedCtrlPts(dim);
    for(unsigned int i = 0; i< dim; i++){
        WeightedCtrlPts[i]=WeightCtrlPts(CtrlPts[i], Weights.GetRow(i));
    }
    return WeightedCtrlPts;
}

/******************************************************************************* */

/******************************************************************************* */
std::vector<iMatrix<double>> UnWeightCtrlPtsSurf(std::vector<iMatrix<double>> &CtrlPtsW){
    int NumPts2= CtrlPtsW.size();
    int NumPts1 = CtrlPtsW[0].GetNumCols();
    int dim = CtrlPtsW[0].GetNumRows()-1;
    
    std::vector<iMatrix<double>> CtrlPts(NumPts2);

    for(unsigned int k=0; k<NumPts2; k++){
        CtrlPts[k]=iMatrix<double>(dim, NumPts1);
        //iMatrix<double> AuxMat (dim,NumPts1);
        
        for(unsigned int j=0; j< NumPts1;j++){
            
            for(unsigned int i=0; i<dim; i++){
                
                if(CtrlPtsW[k](dim,j)==0){
                    
                    //nothing, the declararion automatically sets it to 0
                }
                else{
                   
                    CtrlPts[k](i,j)=CtrlPtsW[k](i,j)/CtrlPtsW[k](dim,j);
                }
            }
        }
       
      

    }
    return CtrlPts;

}

/******************************************************************************* */

/******************************************************************************* */
std::vector<double> CurvePointRational(int n, int p, std::vector<double> &KnotVector, iMatrix<double> &CtrlPtsWeighted, double t){
    int span =FindSpan(KnotVector, t, p, n);
    std::vector<double> Funs = BasisFuns(span,t,p,KnotVector);
    int dim= CtrlPtsWeighted.GetNumRows();
    std::vector<double> Cw(dim, 0.0f);
    for(unsigned int j=0; j<=p; j++){
        //Cw=Cw + Funs[j]*CtrlPtsWeighted[span-p+j]
        for(unsigned int i=0; i<dim; i++){
            Cw[i]=Cw[i]+Funs[j]*CtrlPtsWeighted(i,span-p+j);
        }
    }   
    std::vector<double> C(dim-1);
    double w = Cw[dim-1];

    for(unsigned int i=0; i<dim-1; i++){
        C[i]=Cw[i]/w;
    }
    return C;
}

/******************************************************************************* */

/******************************************************************************* */
int Binomial(int k, int i){
    if (i > k - i) i = k - i;  // C(k, i) = C(k, k-i)
    long long res = 1;
    for (int j = 1; j <= k; j++) {
        res *= k - (i - j);
        res /= j;
    }
    return res;
}

/******************************************************************************* */

/******************************************************************************* */
std::vector<double> RationalBasisFuns(int i, double t, int p, std::vector<double> &KnotVector, std::vector<double> &WeightVector){
    std::vector<double> Funs = BasisFuns(i, t, p, KnotVector);
    int itermax=Funs.size();
    double coef=0.0f;
    int span = FindSpan(KnotVector, t, p, KnotVector.size()-2-p);
            //std::cout<<itermax<<std::endl;
    for(int j=0; j<itermax;j++){
        //std::cout<<"a"<<std::endl;
        coef=coef+Funs[j]*WeightVector[span+j-p];
    }
    for(int j=0; j<itermax;j++){
        //std::cout<<"b"<<std::endl;
        Funs[j]=Funs[j]*WeightVector[span+j-p]/coef;
    }

    return Funs;
}

/******************************************************************************* */

/******************************************************************************* */
iMatrix<double> RatDersBasisFuns(int i, double t, int p, int k,std::vector<double> &KnotVector, std::vector<double> &WeightVector){
    iMatrix<double> R_ders(k+1, p+1);
    iMatrix<double> ders = DerBasisFuns(i, t, p, k ,KnotVector);

    std::vector<double> Wk(k+1, 0.0);
    for (int j = 0; j <= k; ++j){
        for (int i = 0; i < p+1; ++i){
            Wk[j] += WeightVector[i] * ders(j,i);
        }   
    }

    for (int a = 0; a < p+1; ++a) {
    R_ders(0, a)= ders(0,a)*WeightVector[a]/Wk[0]; // (D(0, a) * weights[a]) / Wk[0];
    for (int j = 1; j <= k; ++j) {
            double sum = 0.0;
            for (int m = 1; m <= j; ++m){
                sum += Binomial(j, m) * Wk[m] * R_ders(j - m, a);
                double num = ders(j, a) * WeightVector[a] - sum;
                R_ders(j, a)= num / Wk[0];
            }
                
        }
    }
    
    return R_ders;

}

/******************************************************************************* */

/******************************************************************************* */
iMatrix<double> RatCurveDerivs(iMatrix<double> &Aders, iMatrix<double> &wders, int d){
    int NumCols=Aders.GetNumCols();
    iMatrix<double> CK(d,NumCols);
    //CK[0]=Aders[0]/wders[0]
    std::vector<double> init(NumCols);
    for(unsigned int k=0; k<NumCols; k++){
        init[k]=Aders(0,k)/wders(0,0);
    }   
    CK.SetRow(0, init);
    for(unsigned int k =0; k<=d; k++){
        std::vector<double> v = Aders.GetRow(k);
        for(unsigned int i=1;i<=k;i++){
            int Bin=Binomial(k,i);
            //v = v - Binomial(k, i) * wders[i] * CK[k - i];
            for(unsigned int j=0; j<NumCols; j++){
                v[j]=v[j]-Bin*wders(i,0)*CK(k-i,j);
            }
        }
        //CK[k] = v / wders[0];
        std::vector<double> aux=v;
        for(unsigned int j=0; j<NumCols; j++){
            aux[j]=v[j]/wders(0,0);
        }
        CK.SetRow(k,v);
    }

    return CK;

}

/******************************************************************************* */

/******************************************************************************* */
std::vector<double> SurfacePointRational(int n1, int p1, std::vector<double> &KnotVector1, int n2, int p2,  std::vector<double> &KnotVector2, std::vector<iMatrix<double>> &CtrlPtsWeighted, double t1, double t2){
    int span1=FindSpan(KnotVector1,t1,p1,n1);
    int span2=FindSpan(KnotVector2,t2,p2,n2);
    std::vector<double> Basis1=BasisFuns(span1,t1,p1,KnotVector1);
    std::vector<double> Basis2=BasisFuns(span2,t2,p2,KnotVector2);
    int ind1 = span1-p1;
    int dim = CtrlPtsWeighted[0].GetNumRows();
    std::vector<double> Pt(dim);
    iMatrix<double> auxtemp(p2+1,dim);
    for(unsigned int l=0; l<=p2; l++){

        std::vector<double> temp(dim);
        
        int ind2= span2 - p2+l;
        for(unsigned int k=0;k<=p1;k++){
            double auxBasis1= Basis1[k];
            //std::vector<double> auxCtrlPt = CtrlPts(ind1+k,ind2); 
            for(unsigned int i=0; i<dim;i++){
                //std::cout<<"i="<<i<<",k="<<k<<",l="<<l<<std::endl;
                temp[i]=temp[i]+auxBasis1*CtrlPtsWeighted[ind1+k](i,ind2);
            }
        }
        auxtemp.SetRow(l, temp);


    }

    for(unsigned int l=0; l<=p2;l++){
            double auxBasis2=Basis2[l];
        for(unsigned int i=0; i<dim;i++){
            Pt[i]=Pt[i]+auxBasis2*auxtemp(l,i);
        }
    }

    std::vector<double> C(dim-1);
    double w = Pt[dim-1];

    for(unsigned int i=0; i<dim-1; i++){
        C[i]=Pt[i]/w;
    }
    return C;
}

/******************************************************************************* */

/******************************************************************************* */
iMatrix<std::vector<double>> RatSurfaceDerivs(iMatrix<std::vector<double>> &Aders,iMatrix<double> &wders,int d){

    int dim=Aders(0,0).size();
    iMatrix<std::vector<double>> SKL(d+1,d+1);

    for(unsigned int k=0; k<=d; k++){
        for(unsigned int l=0;l<=d-k;l++){
            //v = Aders [k] [1] ; 
            std::vector<double> v = Aders(k,l);
                                                  
            for(unsigned int j=1; j<=l;j++){
                double binaux=Binomial(1,j);
                double Wdersaux=wders(0,j);
                
                                                  

                for(unsigned int iter=0; iter < dim; iter++){
                    
                    //v = v - Bin[1] [j]*wders[O] [j]*SKL[k] [1-j]; 
                    v[iter]=v[iter]-binaux*Wdersaux*SKL(k,l-j)[iter];

                }
                
                                                 
            }

            for(unsigned int i=1; i<=k; i++){
                double binaux=Binomial(k,i);
                double Wdersaux=wders(i,0);
                
                for(unsigned int iter=0; iter < dim; iter++){
                    //v = v - Bin[k] [i]*wders[i] [O]*SKL[k-i] [1]; 
                    v[iter]=v[iter]-binaux*Wdersaux*SKL(k-i,l)[iter];
                }
                                                    
                //v2 = 0.0;
                std::vector<double> v2(Aders(0,0).size());
                for(unsigned int j=1; j<=l;j++){
                    binaux=Binomial(l,j);
                    Wdersaux=wders(i,j);
                
                    for(unsigned int iter=0; iter < dim; iter++){
                        //v2 = v2 + Bin[l] [j]*wders[i] [j]*SKL[k-i] [l-j]; 
                        v2[iter]=v2[iter]+binaux*Wdersaux*SKL(k-i,l-j)[iter];
                    }
                                                 
                }
                binaux=Binomial(k,i);
                for(unsigned int iter=0; iter < dim; iter++){
                    //v = v - Bin[k] [i]*v2;
                    v[iter]=v[iter]-binaux*v2[iter];
                }
  
                
            }
            //SKL[k] [l] = v/wders[O] [0]; 
            std::vector<double> auxvec(dim);
            for(unsigned int iter=0; iter < dim; iter++){        
                auxvec[iter]=v[iter]/wders(0,0);
            }
                
            SKL(k,l)=auxvec;
        }
    }

    return SKL;
}

/************************************************************************************
*
* Fundamental Geometric Algorithms and support functions
*
************************************************************************************/
int findInsertionIndex(const std::vector<double>& v, double e){
    auto it = std::upper_bound(v.begin(), v.end(), e);
    return it - v.begin()-1;
}

/******************************************************************************* */

/******************************************************************************* */

int countOccurrences(const std::vector<double>& v, double e) {
    double eps = 1e-6f;
    int left = 0;
    int right = v.size() - 1;
    int index = -1;

    //binary search
    while (left <= right) {
        int mid = left + (right - left) / 2;
        if (std::fabs(v[mid] - e) < eps) {
            index = mid;
            right = mid - 1; 
        } else if (v[mid] < e) {
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }

    if (index == -1) return 0; // Not found case

    //linear count
    int count = 0;
    while (index + count < v.size() && std::fabs(v[index + count] - e) < eps) {
        count++;
    }

    return count;
}

/******************************************************************************* */

/******************************************************************************* */
void FindSpanMult(std::vector<double> &KnotVector, double t, int p, int n,int &k, int &s){
    double eps = 1e-6f;
    k=FindSpan(KnotVector,t,p,n);
    s=0;
    for (int i = k; i >= 0 && std::fabs(KnotVector[i] - t) < eps; --i) s++;
    for (int i = k + 1; i <= n + p + 1 && std::fabs(KnotVector[i] - t) < eps; ++i) s++;
}

/******************************************************************************* */

/******************************************************************************* */
iMatrix<double> CurveKnotsIns(int np, int p, std::vector<double> &KnotVector,iMatrix<double> &CtrlPtsW, double u, int index, int s, int r){
    //int mp = np+p+1;
    int dim=CtrlPtsW.GetNumRows();
    int L;
    double alpha;
 
    
    //New CtrlPts
    iMatrix<double> NewCtrlPtsW(dim,CtrlPtsW.GetNumCols()+r);

    for(unsigned int i=0; i<=index-p; i++){
        NewCtrlPtsW.SetCol(i,CtrlPtsW.GetCol(i));
    }

    for(unsigned int i=index-s;i<=np; i++){
        NewCtrlPtsW.SetCol(i+r,CtrlPtsW.GetCol(i)); 
    }
            
    iMatrix<double> AuxCtrlPts(dim,p+1);

    for(unsigned int i=0; i<=p-s;i++){
        AuxCtrlPts.SetCol(i,CtrlPtsW.GetCol(index-p+i));
    }

    for(unsigned int j=1; j<=r;j++){
        
        L=index-p+j;
        
        for(int i=0; i<=p-j-s;i++){
            alpha= (u-KnotVector[L+i])/(KnotVector[i+index+1]-KnotVector[L+i]);
            
            for(unsigned int iter=0; iter<dim; iter++){
                
                AuxCtrlPts(iter,i)=alpha*AuxCtrlPts(iter,i+1)+(1.0f-alpha)*AuxCtrlPts(iter,i);
            }
        }
        
        NewCtrlPtsW.SetCol(L,AuxCtrlPts.GetCol(0));
        
        NewCtrlPtsW.SetCol(index+r-j-s,AuxCtrlPts.GetCol(p-j-s));
        
    }
    for(unsigned int i=L+1;i<index-s; i++){
        
        NewCtrlPtsW.SetCol(i,AuxCtrlPts.GetCol(i-L));
    }

    //Update KnotVector
    KnotVector.insert(KnotVector.begin()+index+1, r, u);

    return NewCtrlPtsW;
}

/******************************************************************************* */

/******************************************************************************* */
std::vector<double> CurvePntByCornerCut(int n, int p, std::vector<double> &KnotVector, iMatrix<double> &CtrlPtsW, double u){
    int dim = CtrlPtsW.GetNumRows()-1;
    std::vector<double> aux(dim+1);
    
    if(u==KnotVector[0]){
        std::vector<double> aux(dim);
        for(unsigned int i=0; i< dim; i++){
            aux[i]=CtrlPtsW(i,0)/CtrlPtsW(dim,0);
        }
        return aux;
    }
    if(u==KnotVector[n+p+1]){
        std::vector<double> aux(dim);
        for(unsigned int i=0; i< dim; i++){
            aux[i]=CtrlPtsW(i,n)/CtrlPtsW(dim,n);
        }
        return aux;
    }
    int k,s;
    FindSpanMult(KnotVector, u, p, n, k, s);
    int r = p-s;
    //std::cout<< "r: " << r << std::endl; 
    iMatrix<double> AuxCtrlPts(dim+1,r+1);
    for(unsigned int i=0; i<=r; i++){
        AuxCtrlPts.SetCol(i,CtrlPtsW.GetCol(k-p+i));
    }

    
    
    for(unsigned int j=1; j<=r; j++){
        for(unsigned int i=0; i<=r-j; i++){
            double alpha= (u-KnotVector[k-p+j+i])/(KnotVector[i+k+1]-KnotVector[k-p+j+i]);
            

            for(unsigned int iter=0; iter< dim+1; iter++){
                aux[iter]=alpha*AuxCtrlPts(iter,i+1)+(1.0f-alpha)*AuxCtrlPts(iter,i);
                
            }
            AuxCtrlPts.SetCol(i,aux);
            //AuxCtrlPts.PrintMatrix();
        }
    }
    std::vector<double> sol(dim);
    
    for(unsigned int i=0; i<dim; i++){
        sol[i]=AuxCtrlPts(i,0)/AuxCtrlPts(dim,0);
    }
    return sol;
}

/******************************************************************************* */

/******************************************************************************* */
std::vector<iMatrix<double>> SurfaceKnotIns(int np1, int p1, std::vector<double> &KnotVector1, int np2, int p2, std::vector<double> &KnotVector2,std::vector<iMatrix<double>> &CtrlPtsW, bool dir,
                                            double uv, int index, int r, int s){
    if(dir==true){//KnotVector1 will be updated
        int dim=CtrlPtsW[0].GetNumRows();
        iMatrix<double> Aux(dim,p1-s+1);
        int numCtrlPts1=CtrlPtsW.size()+1;
        int L;
        iMatrix<double> alpha(p1-s+1,r+1);
        std::vector<iMatrix<double>> NewCtrlPts(numCtrlPts1);
        for(unsigned int i=0; i<numCtrlPts1; i++){
            NewCtrlPts[i]=iMatrix<double>(dim, CtrlPtsW[0].GetNumCols());
        }

        for(unsigned int j=1; j<=p1-j-s;j++){
            L=index-p1+j;
            for(unsigned int i=0;i<=p1-j-s;i++){
                alpha(i,j)= (uv-KnotVector1[L+i])/(KnotVector1[i+index+1]-KnotVector1[L+i]);
            }
        }

        for(unsigned int row=0; row<= np2; row++){
            for(unsigned int i=0; i<=index-p1; i++){
                NewCtrlPts[i].SetCol(row,CtrlPtsW[i].GetCol(row));
            }
            for(unsigned int i=index-s; i<=np1; i++){
                NewCtrlPts[i+r].SetCol(row,CtrlPtsW[i].GetCol(row));
            }
            for(unsigned int i=0; i<=p1-s;i++){
                Aux.SetCol(i,CtrlPtsW[index-p1+i].GetCol(row));
            }
        

            for(unsigned int j=1; j<=r; j++){
                L=index-p1+j;
                for(unsigned int i=0; i<= p1-j-s; i++){
                    std::vector<double> aux(dim);
                    for(unsigned int iter=0; iter<dim; iter++){
                        aux[iter]=  alpha(i,j)*Aux(iter,i+1)+(1.0f-alpha(i,j))*Aux(iter,i);
                    }
                    Aux.SetCol(i,aux);
                }
                NewCtrlPts[L].SetCol(row,Aux.GetCol(0));
                NewCtrlPts[index+r-j-s].SetCol(row,Aux.GetCol(p1-j-s));
            }
            for(unsigned int i=L+1;i<index-s;i++){
                NewCtrlPts[i].SetCol(row,Aux.GetCol(i-L));
            }
        }
        KnotVector1.insert(KnotVector1.begin()+index+1, r, uv);
        return NewCtrlPts;
    }
    else{//KnotVector2 will be updated
    int dim = CtrlPtsW[0].GetNumRows();
    iMatrix<double> Aux(dim, p2 - s + 1);
    int numCtrlPts2 = CtrlPtsW[0].GetNumCols() + 1; 
    int L;
    iMatrix<double> alpha(p2 - s + 1, r + 1);
    std::vector<iMatrix<double>> NewCtrlPts(CtrlPtsW.size());

    for (unsigned int i = 0; i < CtrlPtsW.size(); i++) {
        NewCtrlPts[i] = iMatrix<double>(dim, numCtrlPts2);
    }

    for (unsigned int j = 1; j <= p2 - j - s; j++) {
        L = index - p2 + j;
        for (unsigned int i = 0; i <= p2 - j - s; i++) {
            alpha(i, j)= (uv - KnotVector2[L + i]) / (KnotVector2[i + index + 1] - KnotVector2[L + i]);
        }
    }

    for (unsigned int col = 0; col <= np1; col++) {
        for (unsigned int i = 0; i <= index - p2; i++) {
            NewCtrlPts[col].SetCol(i, CtrlPtsW[col].GetCol(i));
        }
        for (unsigned int i = index - s; i <= np2; i++) {
            NewCtrlPts[col].SetCol(i + r, CtrlPtsW[col].GetCol(i));
        }

        for (unsigned int i = 0; i <= p2 - s; i++) {
            Aux.SetCol(i, CtrlPtsW[col].GetCol(index - p2 + i));
        }

        for (unsigned int j = 1; j <= r; j++) {
            L = index - p2 + j;
            for (unsigned int i = 0; i <= p2 - j - s; i++) {
                std::vector<double> aux(dim);
                for (unsigned int d = 0; d < dim; d++) {
                    aux[d] = alpha(i, j) * Aux(d, i + 1) +
                             (1.0f - alpha(i, j)) * Aux(d, i);
                }
                Aux.SetCol(i, aux);
            }

            NewCtrlPts[col].SetCol(L, Aux.GetCol(0));
            NewCtrlPts[col].SetCol(index + r - j - s, Aux.GetCol(p2 - j - s));
        }

        for (unsigned int i = L + 1; i < index - s; i++) {
            NewCtrlPts[col].SetCol(i, Aux.GetCol(i - L));
        }
    }

    KnotVector2.insert(KnotVector2.begin() + index + 1, r, uv);

    return NewCtrlPts;
        
    }
    
}


/******************************************************************************* */

/******************************************************************************* */
iMatrix<double> RefineKnotVectCurve(int n, int p, std::vector<double> &KnotVector, iMatrix<double> &CtrlPtsW, std::vector<double> &Insertions, int r){
    int m=n+p+1;
    int dim=CtrlPtsW.GetNumRows();
    std::vector<double> NewKnotVector(m+r+2);
    //int numCtrlPts=CtrlPtsW.GetNumCols();
    iMatrix<double> NewCtrlPts(dim,n+r+2);
    int a=FindSpan(KnotVector,Insertions[0],p,n);
    int b=FindSpan(KnotVector,Insertions[r],p,n)+1;

    for(unsigned int j=0;j<=a-p;j++){
        NewCtrlPts.SetCol(j,CtrlPtsW.GetCol(j));
    }   
    for(unsigned int j=b-1;j<=n; j++){
        NewCtrlPts.SetCol(j+r+1,CtrlPtsW.GetCol(j));
    }
    for(unsigned int j=0; j<=a;j++){
        NewKnotVector[j]=KnotVector[j];
    }
    for(unsigned int j=b+p; j<=m; j++){
        NewKnotVector[j+r+1]=KnotVector[j];
    }


    int i=b+p-1;
    int k=b+p+r;
    for(int j=r; j>=0; j--){
        while(Insertions[j]<=KnotVector[i]&&i>a){
            NewCtrlPts.SetCol(k-p-1, CtrlPtsW.GetCol(i-p-1));
            NewKnotVector[k]=KnotVector[i];
            k--;
            i--;
        }
        NewCtrlPts.SetCol(k-p-1,NewCtrlPts.GetCol(k-p));
        for(unsigned int l=1; l<=p;l++){
            int ind = k-p+l;
            double alpha=NewKnotVector[k+l]-Insertions[j];
            if(alpha==0.0){ 
                NewCtrlPts.SetCol(ind-1,NewCtrlPts.GetCol(ind));
            }
            else{
                alpha=alpha/(NewKnotVector[k+l]-KnotVector[i-p+l]);
                std::vector<double> aux(dim);
                for(unsigned int iter=0; iter<dim; iter++){
                    aux[iter]=alpha*NewCtrlPts(iter,ind-1)+(1.0f-alpha)*NewCtrlPts(iter,ind);
                }
                NewCtrlPts.SetCol(ind-1,aux);
            }
            
        }
        NewKnotVector[k]=Insertions[j];
        k--;
    }
    KnotVector=NewKnotVector;
    return NewCtrlPts;
}

/******************************************************************************* */

/******************************************************************************* */
std::vector<iMatrix<double>> RefineKnotVectSurface(int np1, int p1, std::vector<double> &KnotVector1, int np2, int p2, std::vector<double> &KnotVector2,std::vector<iMatrix<double>> &CtrlPtsW, bool dir,
std::vector<double> Insertions, int r){
    int dim=CtrlPtsW[0].GetNumRows();
    if(dir==true){
        //int L;
        int numCtrlPts1=CtrlPtsW.size()+r+1;
        int numCtrlPts2=CtrlPtsW[0].GetNumCols();
        int a=FindSpan(KnotVector1,Insertions[0],p1,np1);
        int b=FindSpan(KnotVector1,Insertions[r],p1,np1)+1;
        int m=np1+p1+1;
        std::vector<double> NewKnotVector(m+r+2);
        std::vector<iMatrix<double>> NewCtrlPts(numCtrlPts1); 
        for(unsigned int i=0; i<numCtrlPts1; i++){
            NewCtrlPts[i]=iMatrix<double>(dim, numCtrlPts2);
        }

        for(unsigned int row=0; row<=m;row++){
            for(unsigned int k=0; k<=a-p1;k++){
                NewCtrlPts[k].SetCol(row,CtrlPtsW[k].GetCol(row));
            }
            for(unsigned int k=b-1; k<=np1; k++){
                NewCtrlPts[k+r+1].SetCol(row,CtrlPtsW[k].GetCol(row));
            }
            

        }
        for(unsigned int j=0; j<=a;j++){
        NewKnotVector[j]=KnotVector1[j];
        }
        for(unsigned int j=b+p1; j<=m; j++){
        NewKnotVector[j+r+1]=KnotVector1[j];
        }
        
        int i=b+p1-1;
        int k=b+p1+r;

        for(int j=r; j>=0; j--){
            while(Insertions[j]<=KnotVector1[i]&&i>a){ 
                for(unsigned int row=0; row<=m;row++){
                    NewCtrlPts[k-p1-1].SetCol(row,CtrlPtsW[i-p1-1].GetCol(row));
                }               
                NewKnotVector[k]=KnotVector1[i];
                k--;
                i--;
            }

            for(unsigned int row=0; row<=m;row++){
                NewCtrlPts[k-p1-1].SetCol(row,NewCtrlPts[k-p1].GetCol(row));
            }

            for(unsigned int l=1; l<=p1; l++){
                int ind=k-p1+l;
                double alpha=NewKnotVector[k+l]-Insertions[j];
                if(alpha==0.0){ 
                    //NewCtrlPts.SetCol(ind-1,NewCtrlPts.GetCol(ind));
                    for(unsigned int row=0; row<=m;row++){
                        NewCtrlPts[ind-1].SetCol(row,NewCtrlPts[ind].GetCol(row));
                    }
                }   
                else{
                    alpha=alpha/(NewKnotVector[k+l]-KnotVector1[i-p1+l]);
                    for(unsigned int row=0; row<=m;row++){
                        std::vector<double> aux(dim);
                        for(unsigned int iter=0; iter<dim; iter++){
                            aux[iter]=alpha*NewCtrlPts[ind-1](iter,row)+(1.0f-alpha)*NewCtrlPts[ind](iter,row);
                        }
                        NewCtrlPts[ind-1].SetCol(row,aux);
                    }
            
                }   
            
            }
        NewKnotVector[k]=Insertions[j];
        k--;
        }

        KnotVector1=NewKnotVector;
        return NewCtrlPts;
    }
    else{
        int numCtrlPts1 = CtrlPtsW.size();
        int numCtrlPts2 = CtrlPtsW[0].GetNumCols() + r + 1;
        int a = FindSpan(KnotVector2, Insertions[0], p2, np2);
        int b = FindSpan(KnotVector2, Insertions[r], p2, np2) + 1;
        int m = np2 + p2 + 1;

        std::vector<double> NewKnotVector(m + r + 2);
        std::vector<iMatrix<double>> NewCtrlPts(numCtrlPts1);

        for (int i = 0; i < numCtrlPts1; i++) {
            NewCtrlPts[i] = iMatrix<double>(dim, numCtrlPts2);
        }

        for (int col = 0; col <= a - p2; col++) {
            for (int i = 0; i < numCtrlPts1; ++i) {
                NewCtrlPts[i].SetCol(col, CtrlPtsW[i].GetCol(col));
            }
        }

        for (int col = b - 1; col <= np2; col++) {
            for (int i = 0; i < numCtrlPts1; i++) {
                NewCtrlPts[i].SetCol(col + r + 1, CtrlPtsW[i].GetCol(col));
            }
        }

        for (int j = 0; j <= a; j++) {
            NewKnotVector[j] = KnotVector2[j];
        }

        for (int j = b + p2; j <= m; j++) {
            NewKnotVector[j + r + 1] = KnotVector2[j];
        }

        int i = b + p2 - 1;
        int k = b + p2 + r;

        for (int j = r; j >= 0; j--) {
            while (Insertions[j] <= KnotVector2[i] && i > a) {
                for (int row = 0; row < numCtrlPts1; ++row) {
                    NewCtrlPts[row].SetCol(k - p2 - 1, CtrlPtsW[row].GetCol(i - p2 - 1));
                }
                NewKnotVector[k] = KnotVector2[i];
                k--;
                i--;
            }

            for (int row = 0; row < numCtrlPts1; row++) {
                NewCtrlPts[row].SetCol(k - p2 - 1, NewCtrlPts[row].GetCol(k - p2));
            }

            for (int l = 1; l <= p2; l++) {
                int ind = k - p2 + l;
                double alpha = NewKnotVector[k + l] - Insertions[j];
                if (alpha == 0.0f) {
                    for (int row = 0; row < numCtrlPts1; row++) {
                        NewCtrlPts[row].SetCol(ind - 1, NewCtrlPts[row].GetCol(ind));
                    }
                } else {
                    alpha = alpha / (NewKnotVector[k + l] - KnotVector2[i - p2 + l]);
                    for (int row = 0; row < numCtrlPts1; row++) {
                        std::vector<double> aux(dim);
                        for (int iter = 0; iter < dim; iter++) {
                            aux[iter] = alpha * NewCtrlPts[row](iter, ind - 1) + (1.0f - alpha) * NewCtrlPts[row](iter, ind);
                        }
                        NewCtrlPts[row].SetCol(ind - 1, aux);
                    }
                }
            }

            NewKnotVector[k] = Insertions[j];
            k--;
        }

        KnotVector2 = NewKnotVector;
        return NewCtrlPts;

    }

}

/************************************************************************************
*
* Functions for the Assembly of the Matrix
*
************************************************************************************/
void legendre_pol(int n, Eigen::VectorXd &nodes, Eigen::VectorXd &weights) {
    
    Eigen::VectorXd beta(n-1);
    for (int i = 0; i < n-1; ++i) {
        double k = i + 1;
        beta(i) = k / std::sqrt(4.0 * k * k - 1);
    }

    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(n, n);
    for (int i = 0; i < n-1; ++i) {
        J(i, i+1) = beta(i);
        J(i+1, i) = beta(i);
    }
 
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(J);
    nodes = solver.eigenvalues();

    Eigen::MatrixXd V = solver.eigenvectors();
    weights.resize(n);
    for (int i = 0; i < n; ++i) {
        weights(i) = 2.0 * V(0, i) * V(0, i);
    }
}


/******************************************************************************* */

/******************************************************************************* */
std::vector<double> removeDuplicates(const std::vector<double>& input) {
    std::vector<double> result = input;
    auto last = std::unique(result.begin(), result.end());
    result.erase(last, result.end());
    return result;
}

/******************************************************************************* */

/******************************************************************************* */
std::vector<double> createSubIntervals(const double nElements){
    std::vector<double> sol(nElements+1);
    double Step=1/nElements;
    for(unsigned int i=0; i<=nElements; i++){
        sol[i]=i*Step;
    }
    return sol;
}

/******************************************************************************* */

/******************************************************************************* */
std::vector<double> removeMatches(const std::vector<double>& v1,  const std::vector<double>& v2   ){
    if (v1.empty()) return {};

    double start = v1.front();
    double step  = v1.size() > 1 ? v1[1] - v1[0] : 1;

    std::vector<char> toRemove(v1.size(), 0); 

    for (size_t i = 0; i < v2.size(); ++i) {
        double val = v2[i];
        if (val < start) continue;

        double diff = val - start;
        double idxReal = diff / step;
        long long idx = llround(idxReal);

        if (std::fabs(idxReal - idx) > 1e-9) continue;

        if (idx >= 0 && static_cast<size_t>(idx) < v1.size()) {
            toRemove[idx] = 1;
        }
    }

    std::vector<double> result;
    result.reserve(v1.size());
    for (size_t i = 0; i < v1.size(); ++i) {
        if (!toRemove[i]) result.push_back(v1[i]);
    }

    return result;
}


/******************************************************************************* */

/******************************************************************************* */
void D1_element_eval(int n, int i, int p, int nEvals, double lower_limit, double upper_limit, std::vector<double> &KnotVector, std::vector<double> &Weights,
                    Eigen::VectorXd &nodes, iMatrix<double> &CtrlPtsW, std::function<double(double)> f,
                    iMatrix<double> &BasisFunsEvals, iMatrix<double> &DerBasisFunsEvals,std::vector<double>& JacobianEvals, std::vector<double>& InverseJacobianEvals, std::vector<double>& funcEvals){
    BasisFunsEvals.Resize(p+1,nEvals);                        
    DerBasisFunsEvals.Resize(p+1,nEvals);
    JacobianEvals.resize(nEvals);
    InverseJacobianEvals.resize(nEvals);
    funcEvals.resize(nEvals);
    double auxeval1= (upper_limit - lower_limit)/2;
    double auxeval2= (upper_limit + lower_limit)/2;

    for(unsigned int k=0; k<nEvals; k++){
        //std::cout << k << std::endl;
        double EvalPoint = auxeval1*nodes(k)+auxeval2;
        //std::cout <<"EvalPoint: " << EvalPoint<< std::endl;
        iMatrix<double> AuxMat1 = RatDersBasisFuns(i,EvalPoint, p, 1, KnotVector, Weights);
        BasisFunsEvals.SetCol(k, AuxMat1.GetRow(0));
        DerBasisFunsEvals.SetCol(k, AuxMat1.GetRow(1));
        //std::cout << "i1" << std::endl;
        iMatrix<double> Aders, wders;
        iMatrix<double> Allders = CurveDerivsAlg1(n, p, KnotVector, CtrlPtsW, EvalPoint, 2);
        //std::cout << "i2" << std::endl;
        int auxInt=Allders.GetNumCols()-1;
        Aders = Allders.GetSubMat(0,1,0, auxInt-1);
        wders = Allders.GetSubMat(0,1,auxInt, auxInt);
        iMatrix<double> AuxMat2 = RatCurveDerivs(Aders, wders, 2);
        //std::cout << "i3" << std::endl;

        std::vector<double> AuxVec= AuxMat2.GetRow(1);
        double AuxVal=0;
        for(unsigned int i=0; i< AuxVec.size(); i++){
            AuxVal+= AuxVec[i]*AuxVec[i];
        }
        AuxVal=sqrt(AuxVal);
        //std::cout << "i4" << std::endl;
        JacobianEvals[k]=AuxVal;
        InverseJacobianEvals[k]=1/AuxVal;
        funcEvals[k]=f(EvalPoint);
    }

}

/******************************************************************************* */

/******************************************************************************* */
iMatrix<double> gauss_legendre_cuadrature_integral_bilinealForm(int p, double lower_limit, double upper_limit, int nEvals, iMatrix<double> &DerBasisFunsEvals, std::vector<double>& InverseJacobianEvals, Eigen::VectorXd &weights){
    iMatrix<double> sol(p+1,p+1);
    double AuxConstant = (upper_limit-lower_limit)/2;

    for(unsigned int i=0; i<=p; i++){
        for(unsigned int j=0; j<=p; j++){
            double AuxVal=0;
            for(unsigned int iter=0; iter<nEvals; iter++){
                AuxVal+=weights[iter]*DerBasisFunsEvals(i, iter)*DerBasisFunsEvals(j,iter)*InverseJacobianEvals[iter];
            }
            AuxVal=AuxConstant*AuxVal;

            sol(i,j)=AuxVal;
            //sol(j,i,AuxVal);
        }
    }
    return sol;
}

/******************************************************************************* */

/******************************************************************************* */
std::vector<double> gauss_legendre_cuadrature_integral_linealForm(int p, double lower_limit, double upper_limit, int nEvals, iMatrix<double> &BasisFunsEvals, std::vector<double>& JacobianEvals, std::vector<double>& funcEvals, Eigen::VectorXd &weights){
    std::vector<double> sol(p+1);
    double AuxConstant = (upper_limit-lower_limit)/2;

    for(unsigned int i=0; i<=p; i++){
        double AuxVal=0;
        for(unsigned int iter=0; iter<nEvals; iter++){
                AuxVal+=weights[iter]*BasisFunsEvals(i, iter)*JacobianEvals[iter]*funcEvals[iter];
            }
        sol[i]=AuxConstant*AuxVal;
    }
    return sol;
}


/******************************************************************************* */

/******************************************************************************* */
void impose_Dirichlet_Condition(int p, int size, Eigen::SparseMatrix<double> &global, Eigen::VectorXd &LinearForm, bool Modify0, bool Modify1, double ValueAt0, double ValueAt1){
    if(Modify0){
        
        global.coeffRef(0, 0) = 1.0;
        LinearForm(0)=ValueAt0;
        for(unsigned int i=1; i<=p; i++){
            global.coeffRef(i, 0) = 0.0;
            LinearForm(i)-=ValueAt0*global.coeffRef(0, i);
            global.coeffRef(0, i) = 0.0;
        }
    }

    if(Modify1){

        global.coeffRef(size, size) = 1.0;
        LinearForm(size)=ValueAt1;
        for(unsigned int i=size-1; i>=size-p; i--){
            global.coeffRef(size, i) = 0.0;
            LinearForm(i)-=ValueAt1*global.coeffRef(i, size);
            global.coeffRef(i, size) = 0.0;
        }
    }
}

/******************************************************************************* */

/******************************************************************************* */
void impose_Newmann_Condition(int size, Eigen::VectorXd &LinearForm, bool Modify0, bool Modify1, double ValueAt0, double ValueAt1){
    if(Modify0){   
        LinearForm(0)-=ValueAt0;
    }
    if(Modify1){
        LinearForm(size)+=ValueAt1;
    }
}

/******************************************************************************* */

/******************************************************************************* */
void impose_Robin_Condition(int p, int size, Eigen::SparseMatrix<double> &global, Eigen::VectorXd &LinearForm, double alpha, double beta, bool Modify0, bool Modify1, double ValueAt0, double ValueAt1){
    
    double BetaInverse=1/beta;
    
    if(Modify0){

        global.coeffRef(0, 0) -= alpha/beta;
        LinearForm(0)-=ValueAt0/beta;
        
    }

    if(Modify1){

        global.coeffRef(size, size) += alpha/beta;
        LinearForm(0)+=ValueAt1/beta;
    }
}

/******************************************************************************* */

/******************************************************************************* */
void writeFunction(Eigen::VectorXd funEvals, std::vector<double> KnotVector, std::vector<double> WeightVector, int n, int p, double lowerLimit, double upperLimit,  std::string name){
    std::ofstream fichero;
    fichero.open(name);
    double Step=(upperLimit-lowerLimit)/(funEvals.size()-p);
    for(unsigned int i=0; i<=funEvals.size()-p; i++){
        double t=lowerLimit + i*Step;
        int span= FindSpan(KnotVector, t, p, n);
        double AuxVal=0;
        std::vector<double> eval= RationalBasisFuns(span, t , p, KnotVector, WeightVector);
        for(int j=0; j<=p; j++){
            AuxVal+=funEvals[span-p+j]*eval[j];
        }
        fichero << AuxVal << " ";
    }
    fichero.close();
}



/******************************************************************************* */

/******************************************************************************* */
void writeFunction(std::function<double(double)> f, double lowerLimit, double upperLimit, double nSteps, iMatrix<double> CtrlPtsW, std::vector<double> KnotVector,
                    int n, int p, std::string name){
    std::ofstream fichero;
    fichero.open(name);
    double Step=(upperLimit-lowerLimit)/nSteps;
    
    for(unsigned int i=0; i<= nSteps; i++){
        //double auxval = CurvePointRational(n, p, KnotVector, CtrlPtsW, lowerLimit + i*Step)[0];
        //std::cout << "point is: "<< auxval << std::endl;
        fichero << f(lowerLimit + i*Step) << " ";
    }

    fichero.close();
}