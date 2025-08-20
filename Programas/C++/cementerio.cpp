/*
*Input: A NurbsVector, the degree of the functions and the number of steps
*Output: iMatrix in which we save each funtion in each of the rows discretized in N_Steps
*/
iMatrix<float> BasisFunctions_2(std::vector<float> &NurbsVector, int p, float N_Steps);

/******************************************************************************* */

/******************************************************************************* */
std::vector<float> BasisFunctions_point_2(std::vector<float> &NurbsVector, int p, float t){
    std::vector<float> Nurb(p+1);

    //1ª iter, p=0
    auto it = std::upper_bound(NurbsVector.begin(), NurbsVector.end(), t); //componente no nula
    int iter =static_cast<int>(std::distance(NurbsVector.begin(), it))-1;
    Nurb[ 0 ]=1.0;


    //std::cout << t << " " << NurbsVector[ static_cast<int>(std::distance(NurbsVector.begin(), it))]<< " "<<  static_cast<int>(std::distance(NurbsVector.begin(), it)) << std::endl;


    for(unsigned int j=1; j<=p;j++){
        std::vector<float> New_Nurb(NurbsVector.size()-j);
         for(unsigned int i=0; i<j;i++){
            float auxLeft = Nurb[i];
            float auxRight = Nurb[i+1];
            //std::cout << "left " << auxLeft << "//right " <<auxRight<< std::endl;
            if(auxLeft!=0){
                auxLeft=((t-NurbsVector[iter])/(NurbsVector[iter+j]-NurbsVector[iter]))*auxLeft;                
            }
            

            if(auxRight!=0){
                auxRight=((NurbsVector[iter+j+1]-t)/(NurbsVector[iter+j+1]-NurbsVector[iter+1]))*auxRight;               
            }
            
            //std::cout << "left " << auxLeft << "//right " <<auxRight<< std::endl;
            New_Nurb[i]=auxLeft+auxRight;
        }
        Nurb=New_Nurb;   
    }

    return Nurb;
}

else if(stoi(argv[1])==4){
    int p=3;
    float N_Steps = 10000;
    std::vector<float> Nurbs={0,0,0,1,2,3,4,4,5,5,5};
    auto t1 = high_resolution_clock::now();
    iMatrix<float> sol(Nurbs.size()-p,N_Steps+1);
    //std::vector<float> Nurb(NurbsVector.size());
    float Step=1/(N_Steps);
    for(unsigned int i = 0; i<N_Steps+1; i++){
        float t=Step*i*5;
        //std::cout << "Iteración iniciada " << i<< std::endl;
        int span = FindSpan(Nurbs, t, p, Nurbs.size()-1-p);
        std::vector<float> eval= BasisFuns_Re(span, t , p, Nurbs);
        //std::cout<<"llega post Basis" << std::endl;
        if(span!=-1){
            for(unsigned int j = span-p; j<=span;j++){
                
                
                    sol.SetElement(j,i,eval[j-span+p]);
                
                
            }
        }
        
    }
    auto t2 = high_resolution_clock::now();
    auto ms_int = duration_cast<milliseconds>(t2 - t1);
    std::cout << ms_int.count() << "ms\n";
    Write_BasisFunctions(sol, "C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/BasisFuncts");
    std::cout<<"fin 4"<<endl;

}

else if(stoi(argv[1])==5){//test Funciones bases por nudos
    std::vector<float> Nurbs={0,0,0,1,2,3,4,4,5,5,5};
    auto t1 = high_resolution_clock::now();
    iMatrix<float> Basis=BasisFunctions(Nurbs,3,10000);
    auto t2 = high_resolution_clock::now();
    auto ms_int = duration_cast<milliseconds>(t2 - t1);
    std::cout << ms_int.count() << "ms\n";

    //Basis.PrintMatrix();
    Write_BasisFunctions(Basis, "C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/BasisFuncts");

    iMatrix<float> Basis_1=BasisFunctions(Nurbs,1,100);
    iMatrix<float> Derivates=BasisFunctions_Derivates(Basis_1,1,Nurbs,2);
    Write_BasisFunctions(Derivates, "C:/Users/carlo/OneDrive/Escritorio/Uni/TFG/DataFiles/Nurbs/BasisFuncts_Derivates");
    std::cout<<"fin 2"<<endl;

}

else if(stoi(argv[1])==6){
    int p=3;
    std::vector<float> Nurbs={0,0,0,1,2,3,4,4,5,5,5};
    std::vector<float> t_values ={0,0.1,0.54,1.35,1.79,2,3,4,4.5,5};
    for(unsigned int iter=0; iter<t_values.size();iter++){
        int i = FindSpan(Nurbs,t_values[iter], p, Nurbs.size()-1-p);
        std::cout << t_values[iter]<<" t: " << i<< " i"<<std::endl;
        std::vector<float> eval= BasisFuns_Re(i, t_values[iter] , p, Nurbs);
        for(float i : eval){
            std::cout<< i << " ";
        }
        std::cout<<std::endl;
        
    }

}




#include <vector>
#include <iostream>
#include <cmath>

// Encuentra el span (intervalo de knots) donde se encuentra el punto xi
int FindSpan(int n, int p, double xi, const std::vector<double>& knot_vector) {
    if (xi >= knot_vector[n]) return n - 1; // Caso especial: último span
    if (xi <= knot_vector[p]) return p;     // Caso especial: primer span

    int low = p;
    int high = n;
    int mid = (low + high) / 2;

    // Búsqueda binaria para encontrar el span
    while (xi < knot_vector[mid] || xi >= knot_vector[mid + 1]) {
        if (xi < knot_vector[mid]) high = mid;
        else low = mid;
        mid = (low + high) / 2;
    }
    return mid;
}

// Calcula las derivadas k-ésimas de las funciones base B-spline
void DersBasisFuns(
    int span, double xi, int p, int k, 
    const std::vector<double>& knot_vector, 
    std::vector<std::vector<double>>& ders
) {
    std::vector<std::vector<double>> ndu(p + 1, std::vector<double>(p + 1));
    ndu[0][0] = 1.0;

    // Almacenamiento para diferencias izquierda y derecha
    std::vector<double> left(p + 1), right(p + 1);

    // Llenar la tabla de funciones base (algoritmo recursivo)
    for (int j = 1; j <= p; ++j) {
        left[j] = xi - knot_vector[span + 1 - j];
        right[j] = knot_vector[span + j] - xi;
        double saved = 0.0;

        for (int r = 0; r < j; ++r) {
            ndu[j][r] = right[r + 1] + left[j - r];
            double temp = ndu[r][j - 1] / ndu[j][r];

            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }

    // Inicializar derivadas (0-ésima derivada = función base)
    for (int j = 0; j <= p; ++j) {
        ders[0][j] = ndu[j][p];
    }

    // Calcular derivadas de orden superior
    for (int r = 0; r <= p; ++r) {
        int s1 = 0, s2 = 1; // Alternar filas en a
        std::vector<std::vector<double>> a(2, std::vector<double>(p + 1));
        a[0][0] = 1.0;

        for (int i = 1; i <= k; ++i) {
            double d = 0.0;
            int rk = r - i, pk = p - i;
            if (r >= i) {
                a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
                d = a[s2][0] * ndu[rk][pk];
            }
            int j1 = (rk >= -1) ? 1 : -rk;
            int j2 = (r - 1 <= pk) ? i - 1 : p - r;

            for (int j = j1; j <= j2; ++j) {
                a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
                d += a[s2][j] * ndu[rk + j][pk];
            }
            if (r <= pk) {
                a[s2][i] = -a[s1][i - 1] / ndu[pk + 1][r];
                d += a[s2][i] * ndu[r][pk];
            }
            ders[i][r] = d;
            std::swap(s1, s2); // Alternar filas
        }
    }

    // Ajustar multiplicadores
    int r = p;
    for (int i = 1; i <= k; ++i) {
        for (int j = 0; j <= p; ++j) {
            ders[i][j] *= r;
        }
        r *= (p - i);
    }
}

// Ejemplo de uso
int main() {
    int p = 2; // Grado del B-spline
    int n = 5; // Número de funciones base - 1
    int k = 2; // Máximo orden de derivada (k <= p)
    double xi = 0.3; // Punto paramétrico

    // Knot vector (ejemplo uniforme)
    std::vector<double> knot_vector = {0, 0, 0, 1, 2, 3, 3, 3};

    int span = FindSpan(n, p, xi, knot_vector);
    std::vector<std::vector<double>> ders(k + 1, std::vector<double>(p + 1));

    DersBasisFuns(span, xi, p, k, knot_vector, ders);

    // Imprimir resultados
    std::cout << "Derivadas de las funciones base en xi = " << xi << ":\n";
    for (int i = 0; i <= k; ++i) {
        std::cout << "Orden " << i << ": ";
        for (int j = 0; j <= p; ++j) {
            std::cout << ders[i][j] << " ";
        }
        std::cout << "\n";
    }

    return 0;
}


int FindSpan(std::vector<float> &NurbsVector, float t, int p, int n){

    if(NurbsVector[n+1]==t){
        return n;
    }

    else{
        int low =p;
        int high=n+1;
        int mid=(low+high)/2;
        while(t<NurbsVector[mid] || t>=NurbsVector[mid+1]){
            if(t<NurbsVector[mid]){
                high=mid;
            }
            else{
                low=mid;
            }
            mid=(low+high)/2;
        }
        return mid;
    }
}