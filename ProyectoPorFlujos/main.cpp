#include<iostream>
#include<iomanip> //esta directiva sirve para usar los manipuladores setw y setprecision
#include<cmath>
#include<string>
#include<fstream> //esta directiva sirve para definir los archivos de flujo de entrada y salida
#include<cstdlib> //esta directiva sirve para usar system("cls"), exit(1) y demas funciones de usuario
#define MAX 10 //esta directiva establece una macro en el codigo fuente, sustituye todas las apariciones de MAX en el programa por el valor de 10
                //es un mecanismo para definir constantes simbolicas
using namespace std;

int Menu(int op1) //retorna la opcion del menu
{
         cout<<endl;
         cout<<"\t\t\t ÉÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍ»\n";
         cout<<"\t\t\t º        Aroximacion polinomial         º\n";
         cout<<"\t\t\t ÈÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍ¼\n";
         cout<<"\t\t\t º     Resolver por:                     º\n";
         cout<<"\t\t\t º     1. Factorizacion  LU              º\n";
         cout<<"\t\t\t º     2. Eliminacion Gaussiana          º\n";
         cout<<"\t\t\t º     3. Salir.                         º\n";
         cout<<"\t\t\t ÈÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍ¼\n";

         cout<<endl<<"\t Ingrese Opcion  ";
         cout<<endl<<"\t ÄÄÄÄÄÄÄÄÄÄÄÄÄÄ : ";
         cin>>op1;
         return op1;


}
int SubMenu(int op2)//retorna la opcion del submenu
{
         cout<<endl;
         cout<<"\t\t\t ÉÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍ»\n";
         cout<<"\t\t\t º        Factorizacion  LU              º\n";
         cout<<"\t\t\t ÈÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍ¼\n";
         cout<<"\t\t\t º     Resolver por:                     º\n";
         cout<<"\t\t\t º     1. Doolittle                      º\n";
         cout<<"\t\t\t º     2. Crout                          º\n";
         cout<<"\t\t\t º     3. Cholesky                       º\n";
         cout<<"\t\t\t º     4. Salir.                         º\n";
         cout<<"\t\t\t ÈÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍ¼\n";

         cout<<endl<<"\t Ingrese Opcion  ";
         cout<<endl<<"\t ÄÄÄÄÄÄÄÄÄÄÄÄÄÄ : ";
         cin>>op2;
         return op2;

}
void titulo()//muestra el titulo por consola
{
     cout<<"\t\t\t ÉÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍ»\n";
     cout<<"\t\t\t º                  AJUSTE DE UNA CURVA POLINOMIAL                º\n";
     cout<<"\t\t\t º               POR EL CRITERIO DE MINIMOS CUADRADOS             º\n";
     cout<<"\t\t\t ÈÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍ¼\n\n";
}
//usa iostream y fstream
void lectura_datos(ifstream& leer_entrada, double *x,double *y,int N)
//en las funciones los flujos se pasan como argumentos como parametros de llamada por referencia (ifstream&)
//la funcion lee los pares de datos que se encuentran en el archivo de entrada
{

    if (!leer_entrada) //verifica la conexion con el archivo de entrada
    {
    cout << "Apertura fallida\n";
    return;
   }


    for (int i=0;i<N;i++)
        leer_entrada>>x[i];//lee los valores del eje-x del archivo de entrada "leer entrada"

    for (int i=0;i<N;i++)
        leer_entrada>>y[i];//lee los valores del eje-y del archivo de entrada "leer entrada"
}
void SumatoriasXin(double *X,double *x,int N,int n)
//la funcion almacena en un vector X las sumatorias xj^(i) de una cantidad N de valores "x" leidos del archivo de entrada para un ajuste de un polinomio de grado n
{
    for (int i=0;i<2*n+1;i++)//2n+1 me referencia las sumatiorias que se necesitan al ajustar los puntos a un polinomio de grado n
    {
        X[i]=0; //inicializacion en cero (condicion para simular una sumatoria)
        for (int j=0;j<N;j++)
            X[i]=X[i]+pow(x[j],i);//almacena en un vector de tamaño 2n+1 las sumatorias de los valores xj^(i)
                                  //donde i varia desde 0 hasta ((2n+1)-1)  (sumatorias necesarias para formar el arreglo matricial de los coeficientes) y j varia desde 0 hasta (N-1)  (numero de puntos leidos desde el archivo de entrada)

    }
}
void SumatoriasYiXin(double *Y,double *x,double *y,int N,int n)
//la funcion almacena en un vector Y las sumatorias Yjxj^(i) de una cantidad N de valores (x,y) leidos del archivo de entrada para un ajuste de un polinomio de grado n
{
    for (int i=0;i<n+1;i++)
    {
        Y[i]=0;//inicializacion en cero (condicion para simular una sumatoria)
        for (int j=0;j<N;j++)
        Y[i]=Y[i]+pow(x[j],i)*y[j];//almacena en un vector de tamaño n+1 las sumatorias de los valores Yjxj^(i)
                                  //donde i varia desde 0 hasta ((n+1)-1)  (sumatorias necesarias formar el vector de terminos independientes) y j varia desde 0 hasta (N-1)  (numero de puntos leidos desde el archivo de entrada)
    }
}
double **arreglo_matricial_SXin(double **A,double *X,int n )
//Dispone las sumatorias xj^(i) almacenadas en el vector X en un arreglo matricial (matriz aumentada)
//retornando una matriz llenada parcialmente con las entradas de la matriz de coeficientes
{
    A=new double*[n+1]; //n+1 corresponde a las filas para una matriz ajustada a un polinomio de grado n(n+1 corresponde al tamaño de la matriz de coeficientes)
    for (int i=0;i<=n;i++)
    {
        A[i]=new double[n+2];//n+2 corresponde a las n+1 columnasde la matriz de coeficientes y mas 1 una columna correspondiente al vector de terminos independientes
        for (int j=0;j<=n;j++)
        {
            A[i][j]=X[i+j];//la potencia "p" a la que esta elevados los xi en la sigma sumatoria coincide con la suma de las coordenadas de la matriz en la cual se ubica dicha sumatoria

        }

    }
    return A;
}
void arreglo_matricialAumentado_SYiXin(double **A,double *Y,int n )
//Dispone las sumatorias Yjxj^(i) almacenadas en el vector Y en la n-esima columna que no ha sido llenada del arreglo matricial (matriz aumentada)
{
    for (int i=0;i<=n;i++)
        A[i][n+1]=Y[i];//almacena las sumatorias Yjxj^(i) en la columna n+1 al recorres la i-filas  desde 0 hasta n
}
double** submatriz(double **matriz, int orden, int i, int j)
//determina la submatriz complementaria del elemento aij
//es la matriz que se obtiene al suprimir la fila i y la columna j
{
	double **subm;//la submatriz tiene una dimension n-1xn-1
                  //la matriz riene una fila y columna menos al omitir la fila i y columna j correspondientes a la coordenada del elemento aij
	int p, q; //indices para la matriz
	int a = 0, b;//indices para la submatriz
	 //a y b se inicializan en cero pues van a contener al primer elemento en recorrido de la matriz luego de omitir
	 //la i-esima fila y la j-esima columna

	//asignacion de memoria dinamica para una matriz
	subm = new double* [orden - 1];
	//subm es un puntero dispuesto como arreglo de tamaño fila
	//que en el ciclo for especifica que para cada puntero en el arreglo en sí apunta a un arreglo de números de tipo double de longitud columna
	for(p = 0; p < orden; p++) {
		if(p==i) continue; //omite la i-esima fila
			subm[a] = new double[orden - 1];

			b = 0;

		for(q = 0; q< orden; q++) {
				if(q==j) continue;//omite la j-esima columna
			subm[a][b] = matriz[p][q];//la matriz subm almacena los elementos restantes que son recorridos en el ciclo for cuando se omite la i-esima fila
                                        //y la j-esima columna
			b++;
		}
		a++;//incrementa el indice de la fila de la submatriz
	}
	return subm;
}

double determinante(double **matriz, int orden)
//calcula el determinante de una matriz A de orden n por el metodo de cofactores desarrollado por la primer columna
//la condicion del calculo del determinante es que la matriz sea cuadrada
{
	if(orden == 1)//caso base
		return matriz[0][0];//el determinante de un numero es el mismo numero

	int i;
	int det = 0;//la variable det se inicializa en cero(condicion para simular la sumatoria)
	for(i = 0; i < orden; i++)
		det += (pow(-1.0,(int)i)) * matriz[i][0] * determinante(submatriz(matriz, orden, i, 0), orden - 1);
		//el determinante se calcula como la suma de los productos de los elementos matriz[i][0] de la  columna cero por
		//sus respectivos adjuntos (pow(-1.0,(int)i))* determinante(submatriz(matriz, orden, i, 0), orden - 1)
		//la expresion (pow(-1.0,(int)i)) obedece a un patron de signos considerando un tablero de ajedrez en donde para cada elemento
		//ubicado en una celda le corresponde un signo mas o un signo menos

	return det;
}
void mostrar_matriz_A(ofstream& escribir_salida,double **A,int n)
//escribe la matriz aumentada A de orden n en un archivo de salida
{
     escribir_salida<<"\n\n A =\n";
    escribir_salida<<setw(5);
//el manipulador setw llama a la funcion miembro width que  establece el ancho de campo del siguiente item que esta en la salida del codigo
    for (int i=0;i<n;i++)
    {
        escribir_salida<<"|";
        for (int j=0;j<=n;j++){
            escribir_salida<<"\t"<<setw(8)<<A[i][j]<<setw(8);
            }escribir_salida<<"|";
 			escribir_salida<<"\n    ";
    }
}
void mostrar_Coef(ofstream& escribir_salida,double **A,int n)
//escribe la matriz A de coeficentes de orden n en un archivo de salida
{
    escribir_salida<<"\n\n A =\n";
    escribir_salida<<setw(5);
    for (int i=0;i<n;i++)
    {
        escribir_salida<<"|";
        for (int j=0;j<n;j++){
            escribir_salida<<"\t"<<setw(8)<<A[i][j]<<setw(8);
            }escribir_salida<<"|";
 			escribir_salida<<"\n    ";
    }
}
double **Transpuesta(double **A, int n)
//retorna la matriz transpuesta de la matriz A
{
    //declaracion dinamica de la matriz aux(transpuesta)
   double**aux ;
   aux=new double* [MAX];
   for (int i=0;i<n;i++)
   {
       aux[i]=new double [n];
       for(int j=0;j<n;j++)
       {
           aux[i][j]=0;
       }
   }
   //almacena  la matriz traspuesta en aux.

   for (int i=0; i<n; i++)
   {
       for (int j=0; j<n; j++)
       {
           aux[j][i]=A[i][j];//reorganiza la matriz ubicando los elementos de la fila i-esima en la columna j-esima
       }
   }
   return aux;
}
bool Determinante_Cero(double **M, int n)
//verifica si el determinande de una matriz M de orden n es igual a cero
//retornando true
{
            if (determinante(M,n) != 0.000)
            {
                return false;
            }
            else
            {
              return true;
            }

}
bool Es_simetrica(double **M, int n)
//verifica si una matriz M de orden n es simetrica usando la funcion transpuesta
//retorna true cuando los elementos de la fila i-esima de la matriz transpuesta son iguales a los elementos de la columna j-esima de la matriz original

//esta es la condicion para utilizar la factorizacion de cholesky
{
    double **Tr=Transpuesta(M,n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (M[i][j] != Tr[i][j])
                return false;
    return true;
}
bool Pivotes_Cero(double **M, int n)
//verifica si algun pivote ubicado en la diagonal de una matriz M de orden n es igual a cero
//retorna true en caso de encontrar un cero en la diagonal

//esta es la condicion para utilizar las tres factorizaciones LU (Doolittle, Crout y Cholesky)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (M[i][i] != 0)
                return false;
    return true;
}
void Doolitle(double **a, int n,double **l, double **u)
//Factorización de A=LU donde L posee la diagonal de unos

{
int i,j,k;
double sum;
	for(i=0;i<n;i++)//para realizar las operaciones es necesario
                    //llenar la diagonal con unos de la matriz L  para i=j
                    //llenar el un triangulo formado sobre la diagonal de unos  en ceros de la matriz L para i<j
                    //llenar el triangulo formado bajo la diagonal de la matriz U para i>j
		for(j=0;j<n;j++)
			if(i>j)
			{
			  u[i][j]=0;
			}
			else if(i==j)
			{
			  l[i][j]=1;
			}
			else
			{
			  l[i][j]=0;
			}

	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			sum=0;
			if(i<=j)
			{
				for(k=0;k<n;k++)
					if(k!=i)//en una primera instancia esta condicion me llena la primera fila de la matriz U cuando i=0
                            //con los elementos de la primer fila de la matriz A, es decir u[0][j]=a[0][j]


					sum=sum+l[i][k]*u[k][j];
				        u[i][j]=a[i][j]-sum;//considera los elementos resultantes al aplicar las operaciones (necesarias) sobre la fila i-esima del elemento aij para el cual dicho elemento se redujo a cero
                                            //llevando la matriz A a su forma escalonada
			}
			else
			{
				for(k=0;k<n;k++)
					if(k!=j)//en una segunda instancia esta condicion me llena la primera columna considerada bajo el primer 1 de la diagonal principal de  la matriz L cuando j=0
                            //con los elementos de la primer columna que estan debajo del primer pivote divididos por este mismo pivote de la matriz A, es decir l[i][0]=a[i][0]/u[0][0]

						sum=sum+l[i][k]*u[k][j];
					l[i][j]=(a[i][j]-sum)/u[j][j];//considera los elementos que se utilizaron como multiplicadores por medio del cual se aplicaron las operaciones (necesarias) sobre la fila i-esima del elemento aij
                                            //respecto al pivote para  el cual dicho elemento se redujo a cero llevando la matriz A a su forma escalonada
			}
		}
	}
}
void Crout (double **a, int n, double **l, double **u)
//Factorización de A=LU donde U posee la diagonal de unos
{
    int k=0;
    double sum=0;
    for(int i=0; i<n; i++){//  Llena las matrices. Llena una columna de la matriz L y luego una fila de la matriz U.
        sum=0;
            for(int k=0; k<n; k++){
                sum=0;
            for(int j=0; j<i; j++){//  Multiplica la fila de la matriz L con la columna de la matriz U.
                sum=sum+l[k][j]*u[j][i];
            }
            l[k][i]=a[k][i]-sum;// Toma el resultado de la multiplicación y le resta a matriz A.
        }//Llena las columnas de la matriz L.
        sum=0;
       for(int k=0; k<n; k++){//  Multiplica la fila de la matriz L con la columna de la matriz U.
            sum=0;
            for(int j=0; j<i; j++){
                sum=sum+l[i][j]*u[j][k];// Toma el resultado de la multiplicación y le resta a matriz A, luego divide
                // entre una posición dada de la matriz L.
            }
            u[i][k]=(a[i][k]-sum)/l[i][i];
        }//Llena las filas de la matriz U.
    }
}
double **Cholesky (double **A, int n)
//Factorizacion de cholesky por columna
//se considera una variante de la eliminacion Gaussiana
{
   //Se crea la matriz L y se inicializa por defecto en ceros.
   //tendria de partida  los ceros del triangilo por encima de la diagonal principal para columnas > filas
   double**L ;
   L=new double* [MAX];
   for (int i=0;i<n;i++)
   {
       L[i]=new double [n];

   }



   //Se hace un recorrido primero por columnas y luego filas.
   for (int j=0; j<n; j++){
       for(int i=0; i<n; i++){

           //Se plantea el proceso para los elementos en la diagonal principal.
           if (i==j)// en una primer instancia se define bajo esta condicion con i=j=0 el primer elemento de la diagonal principal L[0][0]=sqrt(a[0][0])

            {
               double sum1=0;
               for(int k=0;k<i;k++)
               {
                   sum1+=pow(L[i][k], 2);
               }
               L[i][j]=sqrt(A[i][j]-sum1);//almacena los elementos de la diagonal afectados por  las operaciones
                                            //al reducir a cero los elementos que se encuentran bajo al pivote de la diagonal principal
                                            //los pivotes se actulizan cada vez que van a reducir los elementos que se encuentran bajo el pivote anterior
               sum1=0;
           }


           else if (i>j)//en una segunda instancial se define bajo esta condicion con i>j con i y con j=0 los elementos bajo el primer pivote ubicados en la primer columna
                        //L[i][0]=a[i][0]/L[0][0] con i variando de 0 a 2
            {
               double sum2=0;
               for (int k=0; k<j; k++){
                   sum2+=(L[i][k]*L[j][k]);
               }
               L[i][j] =(A[i][j] - sum2)/ L[j][j];//almacena los multiplicadores que se calculan como la raiz de las operaciones que sufre la fila
                                            //al reducir a cero los elementos que se encuentran bajo al pivote de la diagonal principal
               sum2=0;
           }



       }
   }

   //La función retorna la matriz L.
   return L;
}
void pivoteo_parcial(double **A, int n)
//calcula el pivoteo parcial al intercambiar
//las filas del pivote original con el menor pivotante parcial hallado al comparar el numero de mayor valor absoluto
//sobre la columna sobre la cual se encuentra el pivote
{
    for (int i=0;i<n;i++)
        for (int k=i+1;k<n;k++)//recorre los elementos bajo el pivote
            if (abs(A[i][i])<abs(A[k][i]))//verifica si algun elemento bajo el pivote tiene mayor valor absoluto en comparacion con los elementos restantes de la columna sobre la cual se halla el pivote
                for (int j=0;j<=n;j++)
                {
                    double aux=A[i][j];//intercambia la fila del pivote original con la fila del elemento de mayor valor absoluto
                    A[i][j]=A[k][j];
                    A[k][j]=aux;
                }
}
void eliminacion_gaussiana(double **A, int n)
//calcula la eliminacion gaussiana al llervar la matriz A de orden n a su forma escalonada por filas
{
    for (int i=0;i<n-1;i++)            //bucle que realiza la eliminacion gaussiana
        for (int k=i+1;k<n;k++)
            {
                double m=A[k][i]/A[i][i];// multiplicadores para la i-esima columna
                for (int j=0;j<=n;j++)
                    A[k][j]=A[k][j]-m*A[i][j];    //elementos resultantes al operar los elementos bajo el pivote  reduciendolos a cero

            }
}
double **coeficientes(double **A, int n)
//abstrae la matriz de coeficientes de la matriz aumentada
{
    double **a;
    a=new double*[n];
    for ( int i = 0 ; i < n ; i ++)
    {
        a[i]=new double[n];
        for (int j = 0 ; j < n ; j ++)
        {
        a[i][j]=A[i][j] ;
        }
    }
    return a;
}
double *terminos_independientes(double **A, int n)
//abstrae los terminos independientes de la matriz aumentada
{
    double *b=new double[n];
    for (int i=0; i<n; i++)
    {
        b[i] = A[i][n];
    }
 return b;
}
double *Sustitucion_Progresiva (double **A, double *b, int n)
{
   // Define el vector solución de retorno de la función.
   double *y = new double[n];

   // Recorre la matriz por filas de arriba para abajo.
   for (int j=0; j<n; j++){

       // Paso base de la sustitución.
       if (j==0){
           y[j] = (b[j]/A[j][j]); //se inicia la sustitucion progresiva con la primer entrada  de la matriz de coeficientes y
                                  //es decir es de la forma y[0]=b[0]/A[0][0]
       }

       // Paso general de la sustitución.
       else{
           double sum=0;
           for ( int k=0; k<j; k++){
               // Realiza la suma de los productos coeficientes-incógnita de las incógnitas ya calculadas.

               sum+= y[k]*A[j][k];
           }

           // Resta al termino independiente la suma obtenida y la divide por el coeficiente correspondiente que acompaña a la variable en la diagonal.
           //donde la resta simula la trasnposicion de terminos cuando se quiere despejar la variablie que se encuentra en la diagonal
           y[j] = (b[j] - sum)/A[j][j];

       }

   }

   return y;
}
double *Sustitucion_Regresiva (double **M, double *Y, int n)
{
   //  Define el vector solución de retorno de la función.
   double *x = new double[n];

   // Recorre la matriz por filas de abajo para arriba.
   for (int j=n-1; j>=0; j--){

       // Paso base de la sustitución.
       if (j==n-1){
           x[j] = (Y[j]/M[j][j]);//se inicia la sustitucion regrsiva con la ultima entrada  de la matriz de coeficientes
                                  //es decir es de la forma y[0]=b[0]/A[0][0]
       }

       // Paso general de la sustitución.
       else{
           double sum=0;
           // Realiza la suma de los productos coeficientes-incógnita de las incógnitas ya calculadas.
           for ( int k=n-1; k>j; k--){
               sum+= x[k]*M[j][k];
           }
           // Resta al termino independiente la suma obtenida y la divide por el coeficiente correspondiente.
           x[j] = (Y[j] - sum)/ M[j][j];
       }

   }

   return x;
}
double **Multiplicaion_Matriz(double** A,double** B,int Af,int Ac,int Bf,int Bc)
 {
    double **M; //se inicializa la matriz en ceros
                //para el cual en cada elemento aij de la posicion del cero correspondiente se almacenara el producto punto del vector de la fila i-esima
                //por la columna j-esima
    M=new double*[Af];
    for(int i=0;i<Af;i++)
    {
        M[i]=new double[Bc];
        for(int j=0;j<Bc;j++)
        {
            M[i][j]=0;
        }
    }
    //se utiliza 3 bucle for
    //los dos bucles mas externos recorren los elementos que almacenan la resultante M[i][j] del producto acumulado de  un elemento de la fila i,columna k con un elemnto de la fila k, columna j
    //el bucle mas interno recorre en paralelo la columna de la primer matriz con la fila de la segunda matriz
    for(int i=0;i<Af;i++)
    {
        for(int j=0;j<Bc;j++)
        {
            for(int k=0;k<Ac;k++)
            {
                M[i][j]+=A[i][k]*B[k][j];
            }
        }
    }
   return M;
 }
 void mostrar_matriz_L(ofstream& escribir_salida,double **L,int n)
 //escribe la matriz triangular inferior L de orden n en un archivo de salida
 {
    escribir_salida<<"\n\n Matriz triangular inferior ";
	escribir_salida<<"\n\n L =\n";
	escribir_salida<<setw(5);
	for(int i=0;i<n;i++)
	{
	    escribir_salida<<"|";
		for(int j=0;j<n;j++)
        {
           escribir_salida<<"\t"<<setw(8)<<L[i][j]<<setw(8);
		}
            escribir_salida<<"|";
 			escribir_salida<<"\n    ";
	}
 }
 void mostrar_matriz_Lt(ofstream& escribir_salida,double **Lt,int n)
 //escribe la matriz triangular inferior transpuesta Lt de orden n en un archivo de salida
 {
    escribir_salida<<"\n\n Matriz triangular inferior transpuesta ";
	escribir_salida<<"\n\n Lt =\n";
	escribir_salida<<setw(5);
	for(int i=0;i<n;i++)
	{
	    escribir_salida<<"|";
		for(int j=0;j<n;j++)
        {
            escribir_salida<<"\t"<<setw(8)<<Lt[i][j]<<setw(8);
        }
            escribir_salida<<"|";
 			escribir_salida<<"\n    ";
	}
 }
 void mostrar_matriz_U(ofstream& escribir_salida,double **U,int n)
 //escribe la matriz triangular superior U de orden n en un archivo de salida
 {
    escribir_salida<<"\n\n ";
    escribir_salida<<"\n\n Matriz triangular superior ";
	escribir_salida<<"\n\n U =\n";
	escribir_salida<<setw(5);
	for(int i=0;i<n;i++)
	{
	    escribir_salida<<"|";
		for(int j=0;j<n;j++)
        {
          escribir_salida<<"\t"<<setw(8)<<U[i][j]<<setw(8);
        }
            escribir_salida<<"|";
 			escribir_salida<<"\n    ";
	}
 }
  void mostrar_LxU(ofstream& escribir_salida,double **LxU,int n)
 //escribe la matriz triangular superior U de orden n en un archivo de salida
 {
    escribir_salida<<"\n\n ";
    escribir_salida<<"\n\n Matriz produto ";
	escribir_salida<<"\n\n LxU =\n";
	escribir_salida<<setw(5);
	for(int i=0;i<n;i++)
	{
	    escribir_salida<<"|";
		for(int j=0;j<n;j++)
        {
          escribir_salida<<"\t"<<setw(8)<<LxU[i][j]<<setw(8);
        }
            escribir_salida<<"|";
 			escribir_salida<<"\n    ";
	}
 }
 void mostrar_LxLt(ofstream& escribir_salida,double **LxU,int n)
 //escribe la matriz triangular superior U de orden n en un archivo de salida
 {
    escribir_salida<<"\n\n ";
    escribir_salida<<"\n\n Matriz produto ";
	escribir_salida<<"\n\n LxLt =\n";
	escribir_salida<<setw(5);
	for(int i=0;i<n;i++)
	{
	    escribir_salida<<"|";
		for(int j=0;j<n;j++)
        {
          escribir_salida<<"\t"<<setw(8)<<LxU[i][j]<<setw(8);
        }
            escribir_salida<<"|";
 			escribir_salida<<"\n    ";
	}
 }
 void mostrar_solucion_Y(ofstream& escribir_salida,double *Y,int n)
 //escribe la solucion intermedia Y al resolver un siterma Ax=b por factorizacion LU en un archivo de salida
 {
    escribir_salida<<"\n\n La solucion intermedia \n";
	escribir_salida<<" con Yc : \n";

	for(int i=0;i<n;i++)
                escribir_salida<<"\n yc"<<"["<<i+1<<"]="<<"\t"<<setw(6)<<Y[i];
                escribir_salida<<"\n\n";
 }
void mostrar_solucion_X(ofstream& escribir_salida,double *X,int n)
//escribe la solucion final X al resolver un siterma Ax=b por factorizacion LU
 {
    escribir_salida<<"\n\n La solucion final \n";
	escribir_salida<<" con Xc :\n";

	for(int i=0;i<n;i++)
    escribir_salida<<"\n xc"<<"["<<i+1<<"]="<<"\t"<<setw(6)<<X[i];
	escribir_salida<<"\n\n";
 }

double error_ajuste(double *x, double *y, double *a, int N, int n)
//calcula el error del ajuste como la raiz de la desviazion cuadratica (sumatoria de la restar de mi
//y-Esperada (calculada con los parametros hallados)con la y[i] de mis puntos dispersos)
{

    double sum =0;
    for (int i = 0; i < N; i++) {
        double yEsperado = 0;
        for (int j = 0; j <= n; j++)
            yEsperado += a[j]* pow(x[i],j);
        double diff = yEsperado - y[i];
        sum += pow( diff, 2);
    }
    sum=sqrt(sum);
   return sum;
}
double error_cuadratico_medio(double *x, double *y, double *a, int N, int n)
//calcula el error cuadratico medio como la raiz de la desviazion cuadratica media  (sumatoria de la restar de mi
//y-Esperada (calculada con los parametros hallados) con la y[i] de mis puntos dispersos  sobre la cantidad de pares de datos)
{

    double sum =0;
    for (int i = 0; i < N; i++) {
        double yEsperado = 0;
        for (int j = 0; j <= n; j++)
            yEsperado += a[j]* pow(x[i],j);
        double diff = yEsperado - y[i];
        sum += pow( diff, 2);
    }
    sum=sqrt(sum/N);
   return sum;
}
void mostrar_parametros_hallados(ofstream& escribir_salida,double *c, int n, int g)
//escribe los parametros hallados al resolver el sistema Ax=b en un archivo de salida
{
     escribir_salida<<"\n los parametros del polinomio de grado "<<g<< " son:\n\n";
    for (int i=0;i<n;i++)
        escribir_salida<<"c ["<<i<<"] ="<<c[i]<<endl;            // Muestra las variables x^0,x^1,x^2,x^3,...
}
void mostrar_polinomio_ajustado(ofstream& escribir_salida,double *c,int n)
//escribe la ecuacion explicita del polinomio de grado n hallado en un archivo de salida
{
    escribir_salida<<"\n El polinomio que mejor se aproxima al conjunto de datos es:\n\n y=";
    for (int i=0;i<n;i++)
        escribir_salida<<" + ("<<c[i]<<")"<<"x^"<<i;
    escribir_salida<<"\n";
}
void mostrar_error_ajuste(ofstream& escribir_salida,double *x, double *y, double *c,int N, int n)
//escribe la solucion intermedia Y al resolver un siterma Ax=b por factorizacion LU
{
   escribir_salida<<"\n El error del ajuste es:\n";
    escribir_salida<<"e= "<<error_ajuste(x,y,c,N,n);
}
 void mostrar_error_cuadratico_medio(ofstream& escribir_salida,double *x, double *y, double *c, int N, int n)
 //escribe la solucion intermedia Y al resolver un siterma Ax=b por factorizacion LU
{
   escribir_salida<<"\n El error cuadratico medio:\n";
   escribir_salida<<"e_med= "<<error_cuadratico_medio(x,y,c,N,n);
}
int main()
{
 int op1,op2;

 char centinela;
 bool salir=false;
 int n,g,N;



    string nombre_archivoEntrada, nombre_archivoSalida;//variable que almacena el nombre de los archivos de entrada y salida
                                                        //que se quieren conectar al flujo
    ifstream leer_entrada;
    ofstream escribir_salida;

     escribir_salida.setf(ios::fixed);//se utilizan manipuladores para configural el formato de los resultados de salida
                                        //que se escribiran en el archivo de salida
     escribir_salida.setf(ios::showpoint);
     escribir_salida.precision(3);
     titulo();

    cout << "\t\t\t   Se leera una cantidad de pares de datos desde un archivo de \n"
          <<"\t\t\t   entrada y luego estos puntos se ajustaran a un polinomio de  \n"
          <<"\t\t\t   grado n y el resultado se colocara en un archivo de salida.\n\n";



        do
      {


    cout << "\t\t\t Ingrese el nombre del archivo de entrada (sin espacios)\n\n";
    cout <<"\t\t\t Archivo de entrada\n";
     cout<<"\t\t\t ÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ : ";
    cin >> nombre_archivoEntrada;

    cout << "\n\n\n\t\t\t Ingrese el nombre del archivo de salida (sin espacios)\n\n";
    cout <<"\t\t\t Archivo de salida\n";
     cout<<"\t\t\t ÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ : ";
    cin >> nombre_archivoSalida;

    cout << "\n\n";
    cout << "\t\t\t Se procedera a leer la matriz desde el archivo  "
          << nombre_archivoEntrada <<"\n"
          <<"\t\t\t y luego se colocaran los resultado  en el archivo "
          << nombre_archivoSalida << endl;


          leer_entrada.open(nombre_archivoEntrada.c_str( ));//la funcion open conecta el archivo de flujo de entrada
                                                            //la funcion c_str realiza la conversion de objetos
                                                            //convierte un objeto string a uno cstring que es el argumento
                                                            //valido para la funcion open
           leer_entrada.clear();//borra el estado de transmision eof(end of file)al leer el archivo de nuevo
                                //para evitar los fallos al intertar leer el archivo de nuevo
           leer_entrada.seekg(0, ios::beg);//Ubinca la posicion del puntero de lectura al comienzo del archivo
    if (leer_entrada.fail( ))//verifica la apertura fallida del archivo de entrada
       {
          cout << "Apertura del archivo de entrada fallido.\n";
          exit(1);//termina del programa inmediatamente
       }

    escribir_salida.open(nombre_archivoSalida.c_str( ));
    if (escribir_salida.fail( ))
       {
          cout << "Apertura del archivo de salida fallido.\n";
          leer_entrada.close( );

          exit(1);
       }
     cout << "\n\n";
     cout<<"\t\t\t Ingrese la cantidad de pares de datos que que se encuentren en el archivo\n\n";
     cout<<"\t\t\t No pares de datos\n";
     cout<<"\t\t\t ÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ : ";
    cin>>N;
    double *x=new double[N];
    double *y=new double[N];
    lectura_datos(leer_entrada,x,y,N);
    cout << "\n\n";
    cout<<"\t\t\t "<<char(168)<<" Cual es el grado del polinomio sobre el que desea ajustar los datos"<<char(63)<<"\n\n";
    cout<<"\t\t\t Grado del polinomio\n";
    cout<<"\t\t\t ÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ: ";
    cin>>n;
    g=n;

    double *X=new double[2*n+1]; //almacena las sumatorias sigma de  1 ,..., 2n y el 1 corresponde al numero de datos
    SumatoriasXin(X,x,N,n);
    double **A=arreglo_matricial_SXin(A,X,n);
    double *Y=new double[n+1];
    SumatoriasYiXin(Y,x,y,N,n);
    arreglo_matricialAumentado_SYiXin(A,Y,n);
    n=n+1;
    double **k=coeficientes(A,n);
    if(!Determinante_Cero(k,n))
    {
        system("cls");
         op1=Menu(op1);


      switch(op1)
            {
             case 1 :

                        system("cls");
                        op2=SubMenu(op2);


                        switch(op2)
                    {
                    case 1 :
                        {
                             escribir_salida<<"\t\t\t ____________________________________________________________________  \n";
                             escribir_salida<<"\t\t\t|  ________________________________________________________________  |\n";
                             escribir_salida<<"\t\t\t| |                  SOLUCION DE SISTEMAS LINEALES AX=B            | |\n";
                             escribir_salida<<"\t\t\t| |                  POR FACTORIZACION LU DE DOOLITTLE             | |\n";
                             escribir_salida<<"\t\t\t| |________________________________________________________________| |\n";
                             escribir_salida<<"\t\t\t|____________________________________________________________________|\n\n";

                     escribir_salida<<"\n\n Matriz aumentada:\n";
                     mostrar_matriz_A(escribir_salida,A,n);
                     double **a=coeficientes(A,n);
                     double *b=terminos_independientes(A,n);
                     double**L1 ;
                    L1=new double* [MAX];
                    for (int i=0;i<n;i++)
                    {
                    L1[i]=new double [n];
                    }
                    double**U1 ;
                    U1=new double* [MAX];
                    for (int i=0;i<n;i++)
                    {
                    U1[i]=new double [n];
                    }
                     Doolitle(a,n,L1,U1);
                     double **L1xU1=Multiplicaion_Matriz(L1,U1,n,n,n,n);

                      if(!Pivotes_Cero(U1,n))
                    {
                     double* q=Sustitucion_Progresiva(L1,b,n);
                     double* p=Sustitucion_Regresiva(U1,q,n);
                     mostrar_matriz_L(escribir_salida,L1,n);
                     mostrar_matriz_U(escribir_salida,U1,n);
                     mostrar_LxU(escribir_salida,L1xU1,n);
                     mostrar_solucion_Y(escribir_salida,q,n);
                     mostrar_solucion_X(escribir_salida,p,n);
                     mostrar_parametros_hallados(escribir_salida,p,n,g);
                     mostrar_polinomio_ajustado(escribir_salida,p,n);
                     mostrar_error_ajuste(escribir_salida,x,y,p,N,n);
                     mostrar_error_cuadratico_medio(escribir_salida,x,y,p,N,n);
                     escribir_salida.close();
                    }
                    else
                        {
                        cout<<"\n Alguno de los elementos obtenidos sobre la diagonal de la matriz triangular superior \n";
                        cout<<"\n son  ceros :\n";

                        cout<<"\n\n No es posible Resolver el sistema por este metodo de factorizacion:\n";
                        cout<<"\nIntentelo de nuevo!!\n";
                        }
                        }



                    break;

                    case 2 :
                        {
                             escribir_salida<<"\t\t\t ____________________________________________________________________  \n";
                             escribir_salida<<"\t\t\t|  ________________________________________________________________  |\n";
                             escribir_salida<<"\t\t\t| |                  SOLUCION DE SISTEMAS LINEALES AX=B            | |\n";
                             escribir_salida<<"\t\t\t| |                  POR FACTORIZACION LU DE CROUT                 | |\n";
                             escribir_salida<<"\t\t\t| |________________________________________________________________| |\n";
                             escribir_salida<<"\t\t\t|____________________________________________________________________|\n\n";
                    escribir_salida<<"\n\n Matriz aumentada:\n";
                    mostrar_matriz_A(escribir_salida,A,n);
                    double **a=coeficientes(A,n);
                    double *b=terminos_independientes(A,n);
                    double**L2 ;
                    L2=new double* [MAX];
                    for (int i=0;i<n;i++)
                    {
                    L2[i]=new double [n];
                    }
                    double**U2 ;
                    U2=new double* [MAX];
                    for (int i=0;i<n;i++)
                    {
                    U2[i]=new double [n];
                    }
                    Crout(a,n,L2,U2);
                    double **L2xU2=Multiplicaion_Matriz(L2,U2,n,n,n,n);
                    if(!Pivotes_Cero(L2,n))
                    {
                    double* q=Sustitucion_Progresiva(L2,b,n);
                    double* p=Sustitucion_Regresiva(U2,q,n);
                    mostrar_matriz_L(escribir_salida,L2,n);
                    mostrar_matriz_U(escribir_salida,U2,n);
                    mostrar_LxU(escribir_salida,L2xU2,n);
                    mostrar_solucion_Y(escribir_salida,q,n);
                    mostrar_solucion_X(escribir_salida,p,n);
                    mostrar_parametros_hallados(escribir_salida,p,n,g);
                    mostrar_polinomio_ajustado(escribir_salida,p,n);
                    mostrar_error_ajuste(escribir_salida,x,y,p,N,n);
                    mostrar_error_cuadratico_medio(escribir_salida,x,y,p,N,n);
                    escribir_salida.close();
                    }
                    else
                        {
                        cout<<"\n Alguno de los elementos obtenidos sobre la diagonal de la matriz triangular inferior \n";
                        cout<<"\n son  ceros :\n";

                        cout<<"\n\n No es posible Resolver el sistema por este metodo de factorizacion:\n";
                        cout<<"\nIntentelo de nuevo!!\n";
                        }
                        }



                    break;
                    case 3 :
                        {
                             escribir_salida<<"\t\t\t ____________________________________________________________________  \n";
                             escribir_salida<<"\t\t\t|  ________________________________________________________________  |\n";
                             escribir_salida<<"\t\t\t| |                  SOLUCION DE SISTEMAS LINEALES AX=B            | |\n";
                             escribir_salida<<"\t\t\t| |                  POR FACTORIZACION LU DE CHOLESKI              | |\n";
                             escribir_salida<<"\t\t\t| |________________________________________________________________| |\n";
                             escribir_salida<<"\t\t\t|____________________________________________________________________|\n\n";
                             double **a=coeficientes(A,n);
                             double *b=terminos_independientes(A,n);
                            if(Es_simetrica(a,n))
                                {


                     escribir_salida<<"\n\n Matriz aumentada:\n";
                     mostrar_matriz_A(escribir_salida,A,n);
                     double** L=Cholesky(a,n);

                     if(!Pivotes_Cero(L,n))
                    {
                     double** Lt=Transpuesta(L,n);
                     double **LxLt=Multiplicaion_Matriz(L,Lt,n,n,n,n);
                     double* q=Sustitucion_Progresiva(L,b,n);
                     double* p=Sustitucion_Regresiva(Lt,q,n);
                     mostrar_matriz_L(escribir_salida,L,n);
                     mostrar_matriz_Lt(escribir_salida,Lt,n);
                     mostrar_LxLt(escribir_salida,LxLt,n);
                     mostrar_solucion_Y(escribir_salida,q,n);
                     mostrar_solucion_X(escribir_salida,p,n);
                     mostrar_parametros_hallados(escribir_salida,p,n,g);
                     mostrar_polinomio_ajustado(escribir_salida,p,n);
                     mostrar_error_ajuste(escribir_salida,x,y,p,N,n);
                     mostrar_error_cuadratico_medio(escribir_salida,x,y,p,N,n);
                     escribir_salida.close();
                    }
                    else
                    {
                       cout<<"\n Alguno de los elementos obtenidos sobre la diagonal de la matriz triangular inferior \n";
                        cout<<"\n son  ceros :\n";

                        cout<<"\n\n No es posible Resolver el sistema por este metodo de factorizacion:\n";
                        cout<<"\nIntentelo de nuevo!!\n";
                    }
                        }
                        else
                        {
                        cout<<"\nla matriz no es simetrica\n";

                        cout<<"\n\n No es posible Resolver el sistema por este metodo de factorizacion:\n";
                        cout<<"\nIntentelo de nuevo!!\n";
                        }
                                }

                    break;



                    case 4 :
                     salir = true;


                    }



             case 2 :
                 {
                             escribir_salida<<"\t\t\t ____________________________________________________________________  \n";
                             escribir_salida<<"\t\t\t|  ________________________________________________________________  |\n";
                             escribir_salida<<"\t\t\t| |                  SOLUCION DE SISTEMAS LINEALES AX=B            | |\n";
                             escribir_salida<<"\t\t\t| |          POR ELIMINACION GAUSSIANA CON PIVOTEO PARCIAL         | |\n";
                             escribir_salida<<"\t\t\t| |________________________________________________________________| |\n";
                             escribir_salida<<"\t\t\t|____________________________________________________________________|\n\n";

                        escribir_salida<<"\n Matriz aumentada:\n";
                        mostrar_matriz_A(escribir_salida,A,n);
                        pivoteo_parcial(A,n);
                        escribir_salida<<"\nPivotasion parcial de las filas de la Matriz es:\n";
                        mostrar_matriz_A(escribir_salida,A,n);
                        eliminacion_gaussiana(A,n);
                        escribir_salida<<"\n\nMatriz despues de la eliminacion gaussiana:\n";
                        mostrar_matriz_A(escribir_salida,A,n);
                        double **a=coeficientes(A,n);
                        double *b=terminos_independientes(A,n);
                        double *c=Sustitucion_Regresiva(a,b,n);
                        mostrar_parametros_hallados(escribir_salida,c,n,g);
                        mostrar_polinomio_ajustado(escribir_salida,c,n);
                        mostrar_error_ajuste(escribir_salida,x,y,c,N,n);
                        mostrar_error_cuadratico_medio(escribir_salida,x,y,c,N,n);
                        escribir_salida.close();
                 }



                    break;



             case 3 :
                     salir = true;
             }
      system("cls");
       if(!salir)
         {
          cout<<"\n\n"<<char(168)<<" Desea continuar "<<char(63)<<"\n";
          cout<<"\t presione c o C para continuar o cualquier otra tecla para salir \n ";
          cout<<"\t ÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ:\n ";
          cin>>centinela;
         }
        system("cls");


     }
     else
     {
        cout << "\n\n";
        cout<<"\t\t\t La matriz  A de nxn cumple las condiciones para un sistema incompatible que son:\n";
        cout<<"\t\t\t 1. El determinante de A = "<<determinante(k,n)<<"\n";
        cout<<"\t\t\t 2. La matriz A es singular \n";
        cout<<"\t\t\t 3. La matriz no se reduce por filas a la matriz identidad In \n";
        break;
     }
      } while((centinela=='c'||centinela=='C')&&!salir);

      //se asegura al cerrar todos los archivos abiertos antes de finalizar el programa
      leer_entrada.close();
      escribir_salida.close();

    cout << "\n\n";
    cout <<"\t\t\t Fin del programa.\n";

}

