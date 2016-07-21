#include <stdio.h>
#include <stdlib.h>
#include "Matrix.h"
#include "MyRand.h"


int main()
{
    srand(time(0));
    int Dimension = 194;
    /**
        Se genera el número de hormigas como multiplo de la poblacion
    **/
    int NumeroHormigas = Dimension ;

     /***
        Declaracion de los parametros (de cuchareo je je)
        rho --> es el factor de evaporación...
     */
     double rho =0.5;
    Matrix * GrafoCiudades = NewMatrix(Dimension, Dimension);
    Vector * CoorX= NewVector(Dimension);
    Vector * CoorY= NewVector(Dimension);
    Vector *Trayecoria = NewVector(Dimension);

    int MaxIteraciones = Dimension*10;
    /***Generar el grafo */

    CargarCiudad("instances/Qatar",GrafoCiudades, CoorX, CoorY);
    //CargarCiudad("instances/Argentina",GrafoCiudades, CoorX, CoorY);
    GenerarGrafo(GrafoCiudades, CoorX, CoorY);
    /**Resolver por medio de Ant Colony*/

    AntColonySolve(GrafoCiudades,Trayecoria, Dimension, NumeroHormigas, rho, MaxIteraciones  );
  //  PrintVector(CoorX);
   // PrintVector(CoorY);




    /****Imprimir la solución*/


    //PrintSolution(GrafoCiudades, Trayecoria, CoorX, CoorY);

    /**Imprimir La treyectoria..**/
    //printf("Trayectoria:\n");
    //for(int i = 0; i < Trayecoria->Size; i++)
      //  printf("%f ", Trayecoria->data[i]);



    FreeMatrix(GrafoCiudades);
    FreeVector(Trayecoria);
    FreeVector(CoorX);
    FreeVector(CoorY);

    return 0;
}
