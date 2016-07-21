#ifndef ANT_COLONY_H_INCLUDED
#define ANT_COLONY_H_INCLUDED
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"
#define TAMPACKET 100
#define PACKETANT 101
#define PACKETANTBESTLOCAL 102
#define FLAG 103
#define ROOT 0
void GenerarGrafo(Matrix *GrafoCiudades, Vector * X, Vector *Y);
void AntColonySolve(Matrix *GrafoCiudades, Vector * XBest, int Dimension, int NumeroHormigas, double rho, double MinFeromona, double MaxFeromona, double FeromonaInicial, double alpha, double beta, int Maxiteraciones );
void PrintSolution(Matrix *GrafoCiudades, Vector *Trayectoria, Vector *X, Vector *Y);
void Desempaquetar(Matrix *Tau, Matrix *Eta, Matrix *Tabu, Matrix *GrafoCiudades, double *Alfa, double *Beta, int *NumeroHormigas);
void RecorrerCiudades(Matrix *Tau, Matrix *Eta, Matrix *Tabu, int Dimension, double *fbest, Matrix* GrafoCiudades, double alpha, double beta, int NumeroHormigas, Vector *XBest);
void EmpaquetarBest(Vector *XBest, double fbest);
int SizePackageBest(Vector *XBest, double fbest);
void CargarCiudad(char * filename,Matrix *GrafoCiudades, Vector *CoorX, Vector * CoorY);
#endif // ANT_COLONY_H_INCLUDED
