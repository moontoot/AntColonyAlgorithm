#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
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
        Number = rand()%(GrafoCiudades->n+1 );
    }
    return Number;
}
void GenerarGrafo(Matrix *GrafoCiudades, Vector * X, Vector *Y)
{
    //~ for(int i = 0; i < X->Size; i++)
    //~ {
//~ 
        //~ X->data[i]  =  GenerarAleatorio(1.0, 10.0);
        //~ Y->data[i]  =   GenerarAleatorio(1.0, 10.0);
    //~ }
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
void RecorrerCiudades(Matrix *Tau, Matrix *Eta, Matrix *Tabu, int Dimension, double *fbest, Matrix* GrafoCiudades, double alpha, double beta, int NumeroHormigas, Vector *XBest)
{
	
	
 /***Para cada hormiga realizar una estimación de la trayectoria, es decir construir una ruta**/
        for(int i = 0 ; i  < NumeroHormigas; i++)
        {
            Vector * L = NewVector(NumeroHormigas);
            //Dado que no se está dibujando entonces no se aplica hasta el elemento n-2
           for(int j =0; j< Dimension-1; j++)
            {
				Vector * tau = NewVector(Dimension);
				Vector * eta = NewVector(Dimension);
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
				FreeVector(tau);
				FreeVector(eta);
 
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
             /****Efecutar elitismo localmente*/
             if( *fbest > L->data[i])
             {
                 cpyMatrixVector(XBest, Tabu, i );
                 *fbest = L->data[i];
             }
             FreeVector(L);
        }
        
}		
unsigned char * EnrutarPaquete(int *Size, Matrix *Tau
, Matrix *Eta, Matrix *Tabu, Matrix *GrafoCiudades, 
double alpha, double beta, int NumeroHormigas, Vector *Sizes, int IndexProcess )
{
	
	  /***Obtener el tamanio de la poblacion..*/
	  
	  if(IndexProcess == 0)
	  {
		 NumeroHormigas = (int )Sizes->data[IndexProcess];
	  } 
	  else if(IndexProcess < Sizes->Size-1)
	  {
		  NumeroHormigas = (int )Sizes->data[IndexProcess+1] - Sizes->data[IndexProcess] ;
	  }else if(IndexProcess == Sizes->Size-1)
	  {
		  NumeroHormigas = (int) Sizes->data[IndexProcess] - Sizes->data[IndexProcess-1] ;
	  }
  		 int increment=0, packsize=0;
         double product	=0;
         
         MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &increment);
         packsize+=9*increment;

         MPI_Pack_size(1, MPI_DOUBLE, MPI_COMM_WORLD, &increment);
         packsize+=2*increment;
         
         product = Tau->m*Tau->n;
         MPI_Pack_size(product, MPI_DOUBLE, MPI_COMM_WORLD, &increment);
         packsize+=increment;
		 product = Eta->m*Eta->n;
         MPI_Pack_size(product, MPI_DOUBLE, MPI_COMM_WORLD, &increment);
         packsize+=increment;
		 
		 product = Tabu->m*NumeroHormigas;
         MPI_Pack_size(product , MPI_DOUBLE, MPI_COMM_WORLD, &increment);
         
         packsize+=increment;
         product = GrafoCiudades->m*GrafoCiudades->n;
         MPI_Pack_size(product, MPI_DOUBLE, MPI_COMM_WORLD, &increment);
         packsize+=increment;
         
         //Primero enviar al tamanio del paquete....
        
         
         //Reservar la memoria del paquete...
         unsigned char *package = (unsigned char *) malloc(packsize);
         int position=0;
         MPI_Pack(&(Tau->n), 1, MPI_INT, package, packsize, &position, MPI_COMM_WORLD );
         MPI_Pack(&(Tau->m), 1, MPI_INT, package, packsize, &position, MPI_COMM_WORLD );
         MPI_Pack(&(Eta->n), 1, MPI_INT, package, packsize, &position, MPI_COMM_WORLD );
         MPI_Pack(&(Eta->m), 1, MPI_INT, package, packsize, &position, MPI_COMM_WORLD );
        
         MPI_Pack(&(NumeroHormigas), 1, MPI_INT, package, packsize, &position, MPI_COMM_WORLD );
         MPI_Pack(&(Tabu->m), 1, MPI_INT, package, packsize, &position, MPI_COMM_WORLD );
        
         MPI_Pack(&(GrafoCiudades->n), 1, MPI_INT, package, packsize, &position, MPI_COMM_WORLD );
         MPI_Pack(&(GrafoCiudades->m), 1, MPI_INT, package, packsize, &position, MPI_COMM_WORLD );
         MPI_Pack(&NumeroHormigas, 1, MPI_INT, package, packsize, &position, MPI_COMM_WORLD );
         
         MPI_Pack(&alpha, 1, MPI_DOUBLE, package, packsize, &position, MPI_COMM_WORLD );
         MPI_Pack(&beta, 1, MPI_DOUBLE, package, packsize, &position, MPI_COMM_WORLD );
         
         MPI_Pack(Tau->data, Tau->n*Tau->m, MPI_DOUBLE, package, packsize, &position, MPI_COMM_WORLD );
         MPI_Pack(Eta->data, Eta->n*Eta->m, MPI_DOUBLE, package, packsize, &position, MPI_COMM_WORLD );
         int indexsubpoblacion = 0;
         if(IndexProcess  > 0)
         indexsubpoblacion =(int) Sizes->data[IndexProcess-1];
         
         MPI_Pack(&(Tabu->data[indexsubpoblacion*Tabu->m]), NumeroHormigas*Tabu->m, MPI_DOUBLE, package, packsize, &position, MPI_COMM_WORLD );
         
         MPI_Pack(GrafoCiudades->data, GrafoCiudades->n*GrafoCiudades->m, MPI_DOUBLE, package, packsize, &position, MPI_COMM_WORLD );
		 *Size = packsize;
		 return package;
}
void AntColonySolve(Matrix *GrafoCiudades, Vector * XBest, int Dimension, int NumeroHormigas, double rho, double MinFeromona, double MaxFeromona, double FeromonaInicial, double alpha, double beta, int Maxiteraciones )
{


    srand(time(0));
    
    int Rank, NumeroProcesos;
	MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
	MPI_Comm_size(MPI_COMM_WORLD, &NumeroProcesos);
    //INICIALIZACION...
    
    
    /****
        Inicialización de parámetros feromonas
    */

    Matrix *Tau = NewMatrix(Dimension, Dimension);
    Matrix *DTau = NewMatrix(Dimension, Dimension);
    Matrix *NTau = NewMatrix(Dimension, Dimension);
    
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
   
    double ValorMinimo = 1e7;
    double fbest= 1e7;
    int Contador = 0;
    double SumaTotalNTau =0;

    int FlagEnvio = 1;
    int IndexBestAnt=0;
    while(Maxiteraciones-- && Contador < 4)
    {
		
		/**
		 * Enviar la información correspondiente al grafo y a la lista tabu
		 * **/
		int packsizeEnv=0;
		/***Particionar el tamanio de la poblacion*/
        Vector *Sizes= NewVector(NumeroProcesos-1);
        int cont=0;
        for(int i = 0; i < Sizes->Size; i++)
        {
			cont+= NumeroHormigas/(NumeroProcesos-1);
			Sizes->data[i] = cont;
		}
		/**Asignar el resto de hormigas..*/
		Sizes->data[Sizes->Size-1] = cont+(NumeroHormigas-cont);
		
		//Enviar el tamanio del paquete
		for(int i = 1; i < NumeroProcesos; i++)
		{
			char unsigned *package = EnrutarPaquete(&packsizeEnv, Tau, Eta, Tabu, GrafoCiudades, alpha, beta, NumeroHormigas, Sizes, i-1);
			//Flag para que el nodo esclavo continue recibiendo
			MPI_Send(&FlagEnvio,1, MPI_INT, i, FLAG, MPI_COMM_WORLD);
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
		//break;
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
    FlagEnvio = 0;
    //Cerrar todos los nodos esclavos
    for(int i = 1; i < NumeroProcesos; i++)
			MPI_Send(&FlagEnvio,1, MPI_INT, i, FLAG, MPI_COMM_WORLD);

     printf("Distancia mínima encontrada %f \n", fbest );


		FreeMatrix(NTau);
		FreeMatrix(Tau);
		FreeMatrix(DTau);
		FreeMatrix(Tabu);
		FreeMatrix(Eta);

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

void Desempaquetar(Matrix *Tau, Matrix *Eta, Matrix *Tabu, Matrix *GrafoCiudades, double *Alfa, double *Beta, int *NumeroHormigas)
{
		int position=0;
		int increment, packsize=0;
		MPI_Recv(&packsize, 1, MPI_INT, ROOT, TAMPACKET, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
		//Reservar la memoria para almacenar el tamanio del paquete..
		unsigned char *package = (unsigned char *) malloc(packsize);
		
		MPI_Recv(package, packsize, MPI_UNSIGNED_CHAR, ROOT, PACKETANT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
		MPI_Unpack(package, packsize, &position, &(Tau->n), 1, MPI_INT, MPI_COMM_WORLD );
		MPI_Unpack(package, packsize, &position, &(Tau->m), 1, MPI_INT, MPI_COMM_WORLD );
		MPI_Unpack(package, packsize, &position, &(Eta->n), 1, MPI_INT, MPI_COMM_WORLD );
		MPI_Unpack(package, packsize, &position, &(Eta->m), 1, MPI_INT, MPI_COMM_WORLD );
		MPI_Unpack(package, packsize, &position, &(Tabu->n), 1, MPI_INT, MPI_COMM_WORLD );
		MPI_Unpack(package, packsize, &position, &(Tabu->m), 1, MPI_INT, MPI_COMM_WORLD );
		MPI_Unpack(package, packsize, &position, &(GrafoCiudades->n), 1, MPI_INT, MPI_COMM_WORLD );
		MPI_Unpack(package, packsize, &position, &(GrafoCiudades->m), 1, MPI_INT, MPI_COMM_WORLD );
		MPI_Unpack(package, packsize, &position, NumeroHormigas, 1, MPI_INT, MPI_COMM_WORLD );
		
		MPI_Unpack(package, packsize, &position, Alfa, 1, MPI_DOUBLE, MPI_COMM_WORLD );
		MPI_Unpack(package, packsize, &position, Beta, 1, MPI_DOUBLE, MPI_COMM_WORLD );
		
		
		Tau->data = (double *) malloc( Tau->n*Tau->m * sizeof(double));
		Eta->data = (double *) malloc( Eta->n*Eta->m * sizeof(double));
		Tabu->data = (double *) malloc( Tabu->n*Tabu->m * sizeof(double));
		GrafoCiudades->data = (double *) malloc( GrafoCiudades->n*GrafoCiudades->m * sizeof(double));
		
		MPI_Unpack(package, packsize, &position, Tau->data, Tau->m*Tau->n, MPI_DOUBLE, MPI_COMM_WORLD );
		MPI_Unpack(package, packsize, &position, Eta->data, Eta->m*Eta->n, MPI_DOUBLE, MPI_COMM_WORLD );
		MPI_Unpack(package, packsize, &position, Tabu->data, Tabu->m*Tabu->n, MPI_DOUBLE, MPI_COMM_WORLD );
		MPI_Unpack(package, packsize, &position, GrafoCiudades->data, GrafoCiudades->m*GrafoCiudades->n, MPI_DOUBLE, MPI_COMM_WORLD );


					
}
int SizePackageBest(Vector *XBest, double fbest)
{
		int increment=0, packsize=0;
		MPI_Pack_size(1, MPI_DOUBLE, MPI_COMM_WORLD, &increment);
		packsize+=increment;
		MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &increment);
		packsize+=increment;
		MPI_Pack_size(XBest->Size, MPI_DOUBLE, MPI_COMM_WORLD, &increment);
		packsize+=increment;
		return packsize;
}

void EmpaquetarBest(Vector *XBest, double fbest)
{

         double product	=0;
         int packsize = SizePackageBest(XBest, fbest);
         
         //Reservar la memoria del paquete...
         unsigned char *package = (unsigned char *) malloc(packsize);
         int position=0;
         MPI_Pack(&fbest, 1, MPI_DOUBLE, package, packsize, &position, MPI_COMM_WORLD );
         MPI_Pack(&(XBest->Size), 1, MPI_INT, package, packsize, &position, MPI_COMM_WORLD );
         MPI_Pack(XBest->data, XBest->Size, MPI_DOUBLE, package, packsize, &position, MPI_COMM_WORLD );
        
         MPI_Send(package, packsize, MPI_UNSIGNED_CHAR, ROOT,PACKETANTBESTLOCAL, MPI_COMM_WORLD);
         
}
void CargarCiudad(char * filename,Matrix *GrafoCiudades, Vector *CoorX, Vector * CoorY)
{
    FILE *filein = fopen(filename, "r");
    int n;
    fscanf( filein,"%d", &n);
    for(int i = 0; i < n; i++)
    {
        int pos;
        fscanf( filein,"%d %lf %lf",&pos, &(CoorY->data[i]), &(CoorX->data[i]) );
    }

    fclose(filein);
}
