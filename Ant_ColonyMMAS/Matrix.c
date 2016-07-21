#include <stdio.h>
#include <stdlib.h>
//#include <mpi.h>
#include "Matrix.h"
Matrix * NewMatrix(int n, int m)
{
    Matrix *M = (Matrix *) malloc( sizeof(Matrix));
    M->data = ( double *) malloc( n*m*sizeof(double));
    M->n=n;
    M->m=m;
    for(int i = 0; i < n*m; i++) M->data[i] = 0;
    return M;
}
Vector * NewVector(int n)
{
    Vector * V = ( Vector *) malloc( sizeof(Vector));
    V->Size=n;
    V->data = (double *) malloc( n*sizeof(double));
    for(int i=0; i < n; i++) V->data[i] = 0;
    return V;
}
void ResetVector(Vector *X){
    for(int i = 0 ; i < X->Size; i++) X->data[i]=0;
}
void cpyRow(Matrix * A, Matrix *B, int Index    )
{
    for(int i = 0; i < A->m; i++)
    A->data[Index*A->m + i ] = B->data[Index*B->m + i ];
}
void cpyVector(Vector *A, Vector *B)
{
    for(int i = 0; i < A->Size; i++)
    A->data[i] = B->data[i];
    A->Size = B->Size;
}
void cpyVectorMatrix(Vector *A, Matrix *B, int Index)
{
    for(int i = 0; i < A->Size; i++)
    B->data[ Index*B->m + i] = A->data[i];
}
void cpyMatrixVector(Vector *A, Matrix *B, int Index)
{
    for(int i = 0; i < A->Size; i++)
        A->data[i] = B->data[ Index*A->Size + i];

}
void ResetMatrix(Matrix *M){
    for(int i = 0; i < M->m * M->n; i++)
        M->data[i]=0;
}
void PrintMatrix(Matrix * M)
{
    for(int i = 0; i < M->n; i++)
    {
        for(int j = 0; j < M->m; j++)
        {
            printf("%.8f ", M->data[i*M->m + j]);
        }
        printf("\n");
    }
}
void cpyMatrix(Matrix *A, Matrix *B)
{
    for(int i = 0; i < A->n; i++)
    {
        for(int j = 0; j < A->m; j++)
        {
            A->data[i*A->m + j] = B->data[i*B->m + j];
        }
    }
}
void PrintVector(Vector *V)
{
    for(int i = 0; i < V->Size; i++)
    printf("%f\n", V->data[i]);
    printf("\n");
}
void FreeMatrix(Matrix *M)
{
    free(M->data);
    free(M);
}
void FreeVector(Vector *V)
{
    free(V->data);
    free(V);
}

void cpyFilaVector(Matrix * A, Vector * V, int Fila, int Size)
{
    for(int i = 0 ; i < Size ; i++)
    {
        V->data[i] = A->data[ Fila*A->m + i ];
    }
}
void cpyColumnaVector(Matrix * A, Vector * V, int Columna)
{
    for(int i = 0 ; i <A->n ; i++)
    {
        V->data[i] = A->data[ i*A->m + Columna ];
    }
}
void SetConstant(Matrix *A, double k)
{
    for(int i=0; i < A->m*A->n; i++)
        A->data[i] = k;

}
//MPI_Datatype MPI_TypeMatrix(Matrix *M)
//{
//	int blocks[3] ={1,1,M->m * M->n};
//	const int nitems=3;
//    MPI_Datatype types[3]={MPI_INT, MPI_INT, MPI_DOUBLE};
//
//    MPI_Datatype mpi_Matrix;
//
//    MPI_Aint address_n, address_m , address_data, baseaddr;
//    MPI_Get_address( M, &baseaddr);
//    MPI_Get_address( &(M->n), &address_n);
//    MPI_Get_address( &( M->m), &address_m);
//    MPI_Get_address( M->data, &address_data);
//
//    MPI_Aint displacements[3];
//    displacements[0] =  (MPI_Aint)0;
//    displacements[1] =   address_n-baseaddr ;
//    displacements[2] =  address_data - baseaddr;
//
//    MPI_Type_create_struct(nitems, blocks, displacements, types, &mpi_Matrix);
//
//    MPI_Type_commit(&mpi_Matrix);
//    return mpi_Matrix;
//}
//MPI_Datatype MPI_TypeVector(Vector *V)
//{
//	int blocks[2] ={1,V->Size};
//	const int nitems=2;
//    MPI_Datatype types[2]={MPI_INT, MPI_DOUBLE};
//
//    MPI_Datatype mpi_Vector;
//
//    MPI_Aint address_Size, address_data, baseaddr;
//    MPI_Get_address( V, &baseaddr);
//    MPI_Get_address( &(V->Size), &address_Size);
//    MPI_Get_address( V->data, &address_data);
//
//    MPI_Aint displacements[2];
//    displacements[0] =  (MPI_Aint)0;
//    displacements[1] = address_data - baseaddr;
//
//    MPI_Type_create_struct(nitems, blocks, displacements, types, &mpi_Vector);
//
//    MPI_Type_commit(&mpi_Vector);
//    return mpi_Vector;
//}
