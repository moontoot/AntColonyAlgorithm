#include <stdio.h>
#include <stdlib.h>
#include "Matrix.h"
#include "MyRand.h"
#include "Ant_Colony.h"
/***
Se genera el grafo de las ciudades donde cada posición indica la distancia, la traza es cero
*/
int Exist(double Number, Vector *V)
{
    for(int i = 0; i < V->Size; i++)
        if(Number == V->data[i])
        {
            return 1;
        }


    return 0;
}
double GetNumber(Vector *V, Matrix *GrafoCiudades)
{
    int Flag=0;
     double Number=0.0;
    while(Exist(Number, V))
    {
        Number = (int)GenerarAleatorio(0.0, GrafoCiudades->n+1);//rand()%(GrafoCiudades->n+1 );
    }
    return Number;
}
void GenerarGrafo(Matrix *GrafoCiudades, Vector * X, Vector *Y)
{
/*    for(int i = 0; i < X->Size; i++)
    {

        X->data[i]  =  (int)GenerarAleatorio(1.0, 10.0);
        Y->data[i]  =  (int)GenerarAleatorio(1.0, 10.0);
    }*/
        /*X->data[0] = 2;
        X->data[1] = 4;
        X->data[2] = 2;
        X->data[3] = 8;
        X->data[4] = 2;
        X->data[5] = 3;
        X->data[6] = 5;
        X->data[7] = 2;
        X->data[8] = 6;
        X->data[9] = 2;

        Y->data[0] = 3;
        Y->data[1] = 6;
        Y->data[2] = 4;
        Y->data[3] = 7;
        Y->data[4] = 9;
        Y->data[5] = 4;
        Y->data[6] = 7;
        Y->data[7] = 5;
        Y->data[8] = 2;
        Y->data[9] = 1;
*/

    /***
        Generar la matriz de distancias
    */
    for(int i = 0; i < GrafoCiudades->n; i++ )
    {
        for(int j = 0; j < GrafoCiudades->m; j++)
        {
            GrafoCiudades->data[i*GrafoCiudades->m + j] = sqrt( pow(X->data[i] - X->data[j] ,2) + pow( Y->data[i] - Y->data[j],2)  );
        }
    }


/*    for(int i =0; i < GrafoCiudades->n; i++)
    {
        for(int j = 0 ; j < i; j++)
        {
                GrafoCiudades->data[i*GrafoCiudades->m + j]=(int)GenerarAleatorio(1.0, 10.0);
                GrafoCiudades->data[j*GrafoCiudades->m + i]=GrafoCiudades->data[i*GrafoCiudades->m + j];
        }
    }*/
}
void GetPij(Vector * Pij,Vector *tau, double alpha, Vector * eta, double beta, Vector *remain)
{
        double Sumatoria = 0;
        for(int i=0; i  <Pij->Size; i++)
        {
            int index = (int) remain->data[i];

            Pij->data[i] = pow(tau->data[index], alpha) * pow(eta->data[index], beta);

            Sumatoria+=Pij->data[i];
        }
        for(int i=0; i  <Pij->Size; i++)
        {
            Pij->data[i]/=Sumatoria;
        }
}

void InicializarListaTabu(Matrix *Tabu)
{
    /**Se genera un vector de las ciudades */
    Vector *Temporal =NewVector(Tabu->m);
    for(int i =0 ; i < Temporal->Size; i++) Temporal->data[i] = i;

    for(int i = 0 ; i < Tabu->n; i++)
    {
            Tabu->data[i*Tabu->m] = -1;
            int IndexRandom =rand()%Temporal->Size;

            while(Tabu->data[i*Tabu->m] ==-1)
            {
                if(Temporal->data[IndexRandom] == -1)
                {
                    IndexRandom =rand()%Temporal->Size;
                    continue;
                }
                Tabu->data[i*Tabu->m] = Temporal->data[IndexRandom];
                Temporal->data[IndexRandom]=-1;
            }
    }
    FreeVector(Temporal);
}
void GetRemain(Vector * remain, Vector *tabuTemporal, int Dimension)
{


    Vector *Temporal = NewVector( Dimension);
    //Eliminar los elementos asignando un número negativo
//    PrintVector(tabuTemporal);
    for(int i = 0; i < tabuTemporal->Size; i++){
      int Indice = tabuTemporal->data[i];
      Temporal->data[ Indice] =  -1;
      //printf("Indice %d\n", Indice);
    }

    /**
        Asignar las posiciones de las ciudades que
        no existen...
    **/
    int cont=0;
    for(int i = 0; cont < remain->Size ; i++)
    {
        if(Temporal->data[i]!=-1)// continue;
        {
            remain->data[cont]=i;
        cont++;
        }

    }

    FreeVector(Temporal);

}
void AntColonySolve(Matrix *GrafoCiudades, Vector * XBest, int Dimension, int NumeroHormigas, double rho, int Maxiteraciones )
{


    srand(time(0));
    //INICIALIZACION...
    /****
        Inicialización de parámetros feromonas
    */

    Matrix *Tau = NewMatrix(Dimension, Dimension);
    Matrix *DTau = NewMatrix(Dimension, Dimension);
    Matrix *NTau = NewMatrix(Dimension, Dimension);
    double MinFeromona=0.1, MaxFeromona=1000;
    double FeromonaInicial = 0.5, alpha = 1, beta  = 5;
    SetConstant(Tau, FeromonaInicial );
    SetConstant(NTau,FeromonaInicial );
    /***
        La matriz eta es la visibilidad expresada como 1 / d_ii
    */
    Matrix *Eta = NewMatrix(Dimension, Dimension);
    for(int i =0; i < Dimension; i++)
    {
        for(int j =0;  j < Dimension; j++)
        {
            Eta->data[i*Eta->m + j] = 1.0/(GrafoCiudades->data[i*Eta->m + j]+0.000001);
        }
    }

    Matrix *Tabu = NewMatrix(NumeroHormigas, Dimension);
    InicializarListaTabu(Tabu);
     /**
    Obtener la constante Q la cual
    */

    double Q = NumeroHormigas;

    /***
        Obtener la longitud de construida de cada ruta
    */
    Vector * L = NewVector(NumeroHormigas);
    double ValorMinimo = 1e7;
    double fbest= 1e7;
    int Contador = 0;
    double SumaTotalNTau =0;

    Vector * tau = NewVector(Dimension);
    Vector * eta = NewVector(Dimension);
    int IndexBestAnt=0;
    while(Maxiteraciones-- && Contador < 4)
    {

            /***
        Generar la lista tabu con los índices de las ciudades, es decir se asigna
        a cada hormiga el ínidice de la respectiva ciudad, no deben existir valores repetidos..
    */


        /***Para cada hormiga realizar una estimación de la trayectoria, es decir construir una ruta**/
        #pragma omp parallel for
        for(int i = 0 ; i  < NumeroHormigas; i++)
        {
            ResetVector(L);
            //Dado que no se está dibujando entonces no se aplica hasta el elemento n-2
           for(int j =0; j< Dimension-1; j++)
            {
               ResetVector(tau);
               ResetVector(eta);
                Vector * tabuTemporal = NewVector(j+1);
                cpyFilaVector( Tabu ,tabuTemporal, i, j+1);
                int IndiceCiudad = tabuTemporal->data[j];

                cpyFilaVector(Tau, tau, IndiceCiudad, Dimension);

                cpyFilaVector(Eta, eta, IndiceCiudad, Dimension);


                /***
                    Generar un vector de indices que no se encentran en la lista tabu
                */

                Vector *remain = NewVector(Dimension - (j+1));

                /***
                    Obtener las ciudades que no se encuentren en la lista tabu
                **/

                GetRemain(remain, tabuTemporal, Dimension);

                /*Obtener la ruta de la hormiga hasta la ruta n*/
                Vector *Pij = NewVector(remain->Size);
                GetPij(Pij, tau, alpha, eta, beta, remain);
                Vector *Indice = NewVector(1);


                Sample(remain,Pij, Indice);

                Tabu->data[i*Tabu->m+j+1]=Indice->data[0];
                FreeVector(Pij);
                FreeVector(Indice);
                FreeVector(remain);
                FreeVector(tabuTemporal);


            }


             /**Obtener los pesos correspondientes a cada ruta...**/
             for(int k = 0; k < (Tabu->m); k++)
             {
                        if(k< (Tabu->m-1))
                        {
                            int IndexI = (int) Tabu->data[i*Tabu->m+(k)];
                            int IndexJ = (int) Tabu->data[i*Tabu->m+(k+1)];
                            L->data[i] += GrafoCiudades->data[IndexI*GrafoCiudades->m + IndexJ];
                        }
                        else
                        {
                            int IndexI = (int) Tabu->data[i*Tabu->m+(k)];
                            int IndexJ = (int) Tabu->data[i*Tabu->m];
                            L->data[i] += GrafoCiudades->data[IndexI*GrafoCiudades->m + IndexJ];
                        }

             }

            /* for(int k = 0; k < (Tabu->m-1); k++)
             {
                        int IndexI = (int) Tabu->data[i*Tabu->m+(k)];
                        int IndexJ = (int) Tabu->data[i*Tabu->m+(k+1)];

                    DTau->data[IndexI*DTau->m + IndexJ] += (Q / L->data[i]);

             }*/
             /****Efecutar elitismo*/
             if(fbest >L->data[i])
             {
                 cpyMatrixVector(XBest, Tabu, i );
                 fbest = L->data[i];
                 IndexBestAnt=i;
             }
               //printf("fbest %f\n", fbest);
        }
        /**Todas la hormigas terminaron su recorrido*/

            //Sólo aumentar el delta de la mejor ruta....
            ResetMatrix(DTau);
            /***
                Obtener el valor de Delta Tau en base a la mejor evaluación...
            */
                for(int k = 0; k < NumeroHormigas; k++)
                {
                    for(int l = 0; l < Dimension-1; l++)
                    {
                        for(int h = 0; h < Dimension-1; h++)
                        {
                            if(XBest->data[l]==Tabu->data[k*Tabu->m + h] && XBest->data[l+1]==Tabu->data[k*Tabu->m+h+1])
                            {
                                int index1 = (int) XBest->data[l];
                                int index2 = (int) XBest->data[l+1];
                                DTau->data[ index1*DTau->m + index2 ] =DTau->data[ index1*DTau->m + index2 ] +( Q / fbest);
                                if(DTau->data[ index1*DTau->m + index2 ] > MaxFeromona) DTau->data[ index1*DTau->m + index2 ] = MaxFeromona;
                                if(DTau->data[ index1*DTau->m + index2 ]< MinFeromona) DTau->data[ index1*DTau->m + index2 ] = MinFeromona;
                                break;
                            }
                        }
                    }
                }

                /*for(int k = 0; k < NumeroHormigas; k++)
                {
                    for(int l = 0; l < Dimension-1; l++)
                    {
                        for(int h = 0; h < Dimension-1; h++)
                        {

                            int CiudadOrigen = XBest->data[l];
                            int CiudadDestino = XBest->data[l+1];

                            if(l == CiudadOrigen && h == CiudadDestino   )
                            {
                                  DTau->data[ l*DTau->m + h ] =DTau->data[ l*DTau->m + h ] +( Q / fbest);

                                if(DTau->data[ l*DTau->m + h ] > MaxFeromona) DTau->data[ l*DTau->m + h ] = MaxFeromona;
                                if(DTau->data[ l*DTau->m + h ]< MinFeromona) DTau->data[ l*DTau->m + h ] = MinFeromona;
                            }
                        }
                    }
                }*/
                //PrintMatrix(DTau);
                //getchar();

        /***
            Realizar el producto  y suma -----> Tau = rho*Tau +DTau
            Evaporación de la feromona
        */
        for(int l = 0 ; l < Tau->n; l++)
        {
            for(int h = 0; h < Tau->m; h++)
            {
                Tau->data[l*Tau->m + h] = rho*Tau->data[l*Tau->m + h] +DTau->data[l*Tau->m + h];
            }
        }
       /***
          Normalizar el valor de Tau
       ****/
         for(int k = 0; k < Dimension; k++)
         {
                ///Obtener la sumatoria de Tau
                double SumatorioNTau = 0;

                for(int l = 0; l < Dimension; l++) SumatorioNTau+= Tau->data[k*Tau->m + l];
                for(int l = 0; l < Dimension; l++)
                {
                    NTau->data[k*NTau->m + l] = Tau->data[k*Tau->m + l]/SumatorioNTau;
                    ///Se aprovecha este ciclo para obtener la sumatoria total de la matriz
                    SumaTotalNTau+=pow((NTau->data[k*NTau->m + l]-0.5),2);

                }
         }

    /****
        Criterio de paro
    ****/
                double CriterioParo = ( (double) (Dimension*Dimension)/4.0) - SumaTotalNTau;
                SumaTotalNTau=0;
                if(ValorMinimo > CriterioParo)
                {

                    Contador = 0;
                    ValorMinimo = CriterioParo;
                }else{
                    Contador++;
                }

    }

     printf("Distancia mínima encontrada %f \n", fbest );

FreeVector(L);
FreeMatrix(NTau);
FreeMatrix(Tau);
FreeMatrix(DTau);
FreeMatrix(Tabu);
FreeMatrix(Eta);
FreeVector(tau);
FreeVector(eta);
}
void PrintSolution(Matrix *GrafoCiudades, Vector *Trayectoria, Vector *X, Vector *Y)
{
    Matrix *Mapa = NewMatrix(Trayectoria->Size,Trayectoria->Size);
    for(int i = 0; i < Mapa->n; i++)
    {
        for(int j = 0;  j < Mapa->m; j++)
        {
            Mapa->data[i*Mapa->m + j]=-1;
        }
    }

    for(int i = 0; i < X->Size; i++)
        {
            int Col = (int)X->data[i];
            int Row  =Y->Size- (int)Y->data[i];

            Mapa->data[Row*Mapa->m + Col  ] = i;
        }

    puts("");

    for(int i = 0; i < Mapa->n; i++)
    {
        for(int j = 0;  j < Mapa->m; j++)
        {
            if(Mapa->data[i*Mapa->m + j]!= -1)
             printf(ANSI_COLOR_RED "%d ",(int)Mapa->data[i*Mapa->m + j] );
            else
             printf(ANSI_COLOR_CYAN"* ");
        }
        printf("\n");
    }
    FreeMatrix(Mapa);
}
void CargarCiudad(char * filename,Matrix *GrafoCiudades, Vector *CoorX, Vector * CoorY)
{
    FILE *filein = fopen(filename, "r");
    if(!filein) {
        puts("Archivo incorrecto");
        exit(-1);
    }
    int n;
    fscanf( filein,"%d", &n);
    for(int i = 0; i < n; i++)
    {
        int pos;
        fscanf( filein,"%d %lf %lf",&pos, &(CoorX->data[i]), &(CoorY->data[i]) );
        printf("%d %f %f\n", pos,CoorX->data[i],CoorY->data[i]);
    }

    fclose(filein);
}
