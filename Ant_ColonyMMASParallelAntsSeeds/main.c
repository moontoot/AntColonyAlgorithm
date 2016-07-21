#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "Matrix.h"
#include "MyRand.h"
#include "Ant_Colony.h"
#define DIMENSION 38
//~ #define DIMENSION 980
#define INSTAN
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
			Vector *XBest = NewVector(Dimension);
			double fbest=1e16;
			int MaxIteraciones = Dimension*10;
			/***Generar el grafo */
			CargarCiudad("instances/Djibouti",GrafoCiudades, CoorX, CoorY);
			//~ CargarCiudad("instances/Qatar",GrafoCiudades, CoorX, CoorY);
			//~ CargarCiudad("instances/Luxembourg",GrafoCiudades, CoorX, CoorY);
			
			GenerarGrafo(GrafoCiudades, CoorX, CoorY);

			/**Resolver por medio de Ant Colony*/
			int Semilla = 11111;
			//Enviar a cada nodo esclavo la información y la semilla..
			
			///////////////////////////////////////////////////////////////////////////////////
			
			/**
			 * Enviar la información correspondiente al grafo y a la lista tabu
			 * **/
			int packsizeEnv=0;
			
			//AntColonySolve(GrafoCiudades,Trayectoria, Dimension, NumeroHormigas, rho, MinFeromona, MaxFeromona, FeromonaInicial, alpha, beta, MaxIteraciones, Semilla  );
			for(int i = 1; i < NumeroProcesos; i++)
			{
				char unsigned *package = EnrutarPaquete(&packsizeEnv, GrafoCiudades, alpha, beta, NumeroHormigas, rho, MinFeromona, MaxFeromona,FeromonaInicial,MaxIteraciones , i-1);
				//Enviar el tamanio del paquete
				MPI_Send(&packsizeEnv,1, MPI_INT, i, TAMPACKET, MPI_COMM_WORLD);
				//Realizar el envio del paquete...
				MPI_Send(package, packsizeEnv, MPI_UNSIGNED_CHAR, i, PACKETANT, MPI_COMM_WORLD);
				free(package);
			}


			//Obtener el tamanio del paquete correspondiente a la solución de cada nodo...
	    //Todos los best tienen el mismo tamanio
	    int SizePacketRecv = SizePackageBest(XBest, fbest);
	    	
		int SizeProcess;
		MPI_Comm_size(MPI_COMM_WORLD, &SizeProcess);
		
		
		double fbestGlobal =1e9;
		Vector *XBestGlobal = NewVector(XBest->Size);
		/***Recibir todos el mejor camino encontrado que corresponde cada subpoblacion*/
		for(int i =1; i  < SizeProcess; i++)
		{
			int position=0;
			unsigned char *PackageRecv = (unsigned char *) malloc(SizePacketRecv);
			MPI_Recv(PackageRecv, SizePacketRecv, MPI_UNSIGNED_CHAR, i, PACKETANTBESTLOCAL, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			
			MPI_Unpack( PackageRecv , SizePacketRecv, &position, &fbestGlobal, 1, MPI_DOUBLE, MPI_COMM_WORLD );

			MPI_Unpack(PackageRecv, SizePacketRecv, &position, &(XBestGlobal->Size), 1, MPI_INT, MPI_COMM_WORLD );

			MPI_Unpack(PackageRecv, SizePacketRecv, &position, XBestGlobal->data, Dimension, MPI_DOUBLE, MPI_COMM_WORLD );
			if(fbestGlobal < fbest )
			{
				fbest = fbestGlobal ;
				cpyVector(XBest, XBestGlobal);
			
			}
     		 free(PackageRecv);
		}
		FreeVector(XBestGlobal); 
		////////////////////////////////////////////////////////////////
			/****Imprimir el mapa en la terminal...*/
			//PrintSolution(GrafoCiudades, XBest, CoorX, CoorY);

			/**Imprimir La treyectoria..**/
		
			printf("Trayectoria:\n");
			for(int i = 0; i < XBest->Size; i++)
				printf("%f ", XBest->data[i]);
			//Imprimir el resultado en R 
			/**Imprimir el mapa en R*/
			PrintR(XBest, CoorX, CoorY);
			
			printf("\nDistancia minima encontrada de todos los procesos %f\n", fbest);	
			FreeMatrix(GrafoCiudades);
			FreeVector(XBest );
			FreeVector(CoorX);
			FreeVector(CoorY);
		}
		else
		{
			
					int Size = NumeroHormigas/NumeroProcesos;
					double Alfa, Beta, fbest=1e9, rho, MinFeromona, MaxFeromona, FeromonaInicial;
					int NumeroHormigas, MaxIteraciones, Semilla;
					Matrix *GrafoCiudades = (Matrix *) malloc(sizeof(Matrix));
					//Recibir y desempaquetar
					Desempaquetar(GrafoCiudades, &Alfa, &Beta, &NumeroHormigas, &rho, &MinFeromona, &MaxFeromona, &FeromonaInicial, &MaxIteraciones);	
					Vector *XBest  = NewVector(GrafoCiudades->m);
					AntColonySolve(GrafoCiudades,XBest, GrafoCiudades->m, NumeroHormigas, rho, MinFeromona, MaxFeromona, FeromonaInicial, Alfa, Beta, MaxIteraciones, Semilla, &fbest  );
					//Empaquetar y enviar al ROOT
					 EmpaquetarBest(XBest, fbest);
					free(GrafoCiudades);
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
			sprintf(Comando, "echo \" plot( x = %s, y = %s,xlab='Horizonte Norte', ylab='Horizonte este', main=c('Trayectoria de las hormigas'), type='p', col='gray', cex=3 );  points(%s, %s,col='red',pch=16); %s  \" | R --Silent --no-save 2>/dev/null | tail -n 0", Nx, Ny, Nx, Ny,Aristas );
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
