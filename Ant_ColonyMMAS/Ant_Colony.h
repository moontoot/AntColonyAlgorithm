#ifndef ANT_COLONY_H_INCLUDED
#define ANT_COLONY_H_INCLUDED
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"
void GenerarGrafo(Matrix *GrafoCiudades, Vector * X, Vector *Y);
void AntColonySolve(Matrix *GrafoCiudades, Vector * Trayectoria, int Dimension, int NumeroHormigas, double rho, int Maxiteraciones );
void PrintSolution(Matrix *GrafoCiudades, Vector *Trayectoria, Vector *X, Vector *Y);
void CargarCiudad(char * filename,Matrix *GrafoCiudades, Vector *CoorX, Vector * CoorY);
#endif // ANT_COLONY_H_INCLUDED
