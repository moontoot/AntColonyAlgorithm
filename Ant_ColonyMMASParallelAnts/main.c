#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "Matrix.h"
#include "Ant_Colony.h"
#include "MyRand.h"
#define DIMENSION 194
void PrintR(Vector * Trayectoria, Vector *CoorX, Vector *CoorY);
int main(int argc, char **argv)
{

	int Rank, NumeroProcesos;
	/***Inicializar el ambiente de MPI*/
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &NumeroProcesos);
	MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
	int Dimension = DIMENSION;
	int NumeroHormigas = Dimension ;
	if(!Rank)
	{
			srand(time(0));


			 /***
				Declaracion de los parametros (de cuchareo je je)
				rho --> es el factor de evaporación...
			 */
			double rho =0.5;
			double MinFeromona=0.1, MaxFeromona=1000;
			double FeromonaInicial = 0.5, alpha = 1, beta  = 5;
			Matrix * GrafoCiudades = NewMatrix(Dimension, Dimension);
			Vector * CoorX= NewVector(Dimension);
			Vector * CoorY= NewVector(Dimension);
			Vector *Trayectoria = NewVector(Dimension);

			int MaxIteraciones = Dimension*10;
			/***Generar el grafo */
			//~ CargarCiudad("instances/Djibouti",GrafoCiudades, CoorX, CoorY);
			CargarCiudad("instances/Qatar",GrafoCiudades, CoorX, CoorY);
			//~ CargarCiudad("instances/Luxembourg",GrafoCiudades, CoorX, CoorY);
			GenerarGrafo(GrafoCiudades, CoorX, CoorY);

			/**Resolver por medio de Ant Colony*/
			AntColonySolve(GrafoCiudades,Trayectoria, Dimension, NumeroHormigas, rho, MinFeromona, MaxFeromona, FeromonaInicial, alpha, beta, MaxIteraciones  );

			/****Imprimir la solución*/


			//PrintSolution(GrafoCiudades, Trayectoria, CoorX, CoorY);

			/**Imprimir La treyectoria..**/

			printf("Trayectoria:\n");
			for(int i = 0; i < Trayectoria->Size; i++)
				printf("%f ", Trayectoria->data[i]);


			//Imprimir el resultado en R

			PrintR(Trayectoria, CoorX, CoorY);


			FreeMatrix(GrafoCiudades);
			FreeVector(Trayectoria );
			FreeVector(CoorX);
			FreeVector(CoorY);
		}
		else
		{

			int Flag=1;
			while(Flag)
			{
					//Comprobar el flag...

					MPI_Recv(&Flag, 1, MPI_INT, ROOT, FLAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

					if(!Flag) break;
					int Size = NumeroHormigas/NumeroProcesos;
					double Alfa, Beta, fbest=1e9;
					int NumeroHormigas;
					Matrix *Tau = (Matrix *) malloc(sizeof(Matrix)) ,
					*Eta = (Matrix *) malloc(sizeof(Matrix)),
					*Tabu = (Matrix *) malloc(sizeof(Matrix)),
					*GrafoCiudades = (Matrix *) malloc(sizeof(Matrix));
					//Recibir y desempaquetar
					Desempaquetar(Tau, Eta, Tabu, GrafoCiudades, &Alfa, &Beta, &NumeroHormigas);
					Vector *XBest  = NewVector(GrafoCiudades->m);
					RecorrerCiudades(Tau, Eta, Tabu, Dimension, &fbest, GrafoCiudades, Alfa, Beta, NumeroHormigas, XBest );
					//Empaquetar y enviar al ROOT
					 EmpaquetarBest(XBest, fbest);
					 
					free(Eta);
					free(Tabu);
					free(GrafoCiudades);
			}
		}
		MPI_Finalize();

    return 0;
}
void PrintR(Vector * Trayectoria, Vector *CoorX, Vector *CoorY)
{

	char *Comando = (char*) malloc(Trayectoria->Size*400* sizeof(char));
	char *Aristas = (char*) malloc(Trayectoria->Size*400* sizeof(char));
	char *Nx = (char*) malloc(Trayectoria->Size*10* sizeof(char));
	char *Ny = (char*) malloc(Trayectoria->Size*10* sizeof(char));
	char *x0 = (char*) malloc(Trayectoria->Size*10* sizeof(char));
	char *x1 = (char*) malloc(Trayectoria->Size*10* sizeof(char));
	char *y0 = (char*) malloc(Trayectoria->Size*10* sizeof(char));
	char *y1 = (char*) malloc(Trayectoria->Size*10* sizeof(char));



	strcpy(Nx, "c(");
	strcpy(Ny, "c(");

	strcpy(x0, "c(");
	strcpy(x1, "c(");
	strcpy(y0, "c(");
	strcpy(y1, "c(");


	for(int i = 0 ; i < Trayectoria->Size; i++)
	{
		int CiudadOrigen = (int)Trayectoria->data[i];
		int CiudadDestino = (int)Trayectoria->data[i+1];
		char temp1[300];
		char temp2[300];
		char tempx0[300];
		char tempx1[300];
		char tempy0[300];
		char tempy1[300];
		if(i==0)
		{
		  sprintf(temp1," %.2f ",CoorX->data[i] );
		  sprintf(temp2," %.2f ",CoorY->data[i] );
		  if(i+1  < Trayectoria->Size )
		  {
			  sprintf(tempx0," %.2f ",CoorX->data[CiudadOrigen] );
			  sprintf(tempx1," %.2f ",CoorX->data[CiudadDestino] );
			  sprintf(tempy0," %.2f ",CoorY->data[CiudadOrigen] );
			  sprintf(tempy1," %.2f ",CoorY->data[CiudadDestino] );
		 }
		}
		else
		{
		  sprintf(temp1,",%.2f ",CoorX->data[i] );
		  sprintf(temp2,",%.2f ",CoorY->data[i] );
		    if(i+1  < Trayectoria->Size )
		  {
		  sprintf(tempx0,",%.2f ",CoorX->data[CiudadOrigen] );
		  sprintf(tempx1,",%.2f ",CoorX->data[CiudadDestino] );
		  sprintf(tempy0,",%.2f ",CoorY->data[CiudadOrigen] );
		  sprintf(tempy1,",%.2f ",CoorY->data[CiudadDestino] );
		}
		}

		strcat(Nx, temp1);
		strcat(Ny, temp2);
		strcat(x0, tempx0);
		strcat(x1, tempx1);
		strcat(y0, tempy0);
		strcat(y1, tempy1);

	}

	strcat(Nx, ")");
	strcat(Ny, ")");
	strcat(x0, ")");
	strcat(x1, ")");
	strcat(y0, ")");
	strcat(y1, ")");

	sprintf(Aristas, " arrows(x0=%s,x1=%s,y0=%s,y1=%s, col='red'); ", x0,x1,y0,y1  );
	//~ Opcion 1 interfacear..
			sprintf(Comando, "echo \" plot( x = %s, y = %s,xlab='Horizonte Norte', ylab='Horizonte este', main=c('Trayectoria de las hormigas'), type='p', col='gray', cex=3 );  points(%s, %s,col='red',pch=16); %s  \" | R --Silent --no-save | tail -n 0", Nx, Ny, Nx, Ny,Aristas );
			system(Comando);
	//~ Opcion 2 scribir un script en R y ejecutarlo
	//~ sprintf(Comando, "echo \" plot( x = %s, y = %s,xlab='Horizonte Norte', ylab='Horizonte este', main=c('Trayectoria de las hormigas'), type='p', col='gray', cex=3 );  points(%s, %s,col='red',pch=16); %s  \" > scripttemporal.R", Nx, Ny, Nx, Ny,Aristas );
	//~ system(Comando);
	//~ system("Rscript scripttemporal.R");
	free(Comando);
	free(Aristas);
	free(Nx);
	free(Ny);
	free(x0);
	free(x1);
	free(y0);
	free(y1);
}
