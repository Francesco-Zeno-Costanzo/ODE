#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define a 1

double* Force(int N, double *x){
	/*
	funzione che prende le posizioni degli oscillatori
	a tempo fisso e restituisce la forza su tutta la catena
	*/
	double *F;
	F = (double *)malloc(N*sizeof(double));

	F[0] = 0;
	for(int k = 1; k < N - 1; k++){
		F[k] = (-2.0*x[k]+x[k - 1]+x[k + 1])*(1.0+a*(x[k + 1]-x[k - 1]));
	}
	F[N - 1] = 0;
	
	return F;
}

double** integratore(int N, int M, double dt){
	/*
	Integratore con il metodo velocity verlet
	*/
	double **X, **V;
	double *F0, *F1, *temp;
	
	X=(double **)malloc(N*sizeof(double));
	for(int i=0; i < N; i++){
    		X[i] = (double*)calloc(M+1, sizeof(double) );
	}
	
	
	V=(double **)malloc(N*sizeof(double));
	for(int i=0; i < N; i++){
    		V[i] = (double*)calloc(M+1, sizeof(double) );
	}
	
	F0 = (double*)malloc(N*sizeof(double));
	F1 = (double*)malloc(N*sizeof(double));
	temp = (double*)malloc(N*sizeof(double));
	
	//condizioni inziali
	for(int i = 0; i < N; i++){
		X[i][0] = sin(0.08975*i);
		V[i][0] = 0;	
	}
	//ciclo sui tempi
	for(int i=0; i < M; i++){
	
		// per ogni tempo ciclo sulle particelle
		for(int l = 0; l < N; l++){
			temp[l] = X[l][i];			
		}
		F0 = Force(N, temp);

		for(int j = 0; j < N; j++){
			X[j][i + 1] = X[j][i] + V[j][i]*dt + F0[j]*((dt*dt)/2.0);
		}
        	for(int l = 0; l < N; l++){
			temp[l] = X[l][i+1];
		}
		
		F1 = Force(N, temp);

    		for(int j = 0; j < N; j++){
			V[j][i + 1] = V[j][i] + (F0[j] + F1[j])*dt/2.0;
		}
	}
	return X;

}

int main(void){
	int N, M;
	double dt;
	N = 36;
	M = 800000;
	dt = 5971/(M*1.0);
	
	double **X;
	X=(double **)malloc(N*sizeof(double));
	for(int i = 0; i < N; i++){
    		X[i] = (double*)calloc(M+1, sizeof(double) );
	}
	X=integratore(N, M, dt);
	
	FILE *fd;
	fd=fopen("scrivi.txt", "w");
	if( fd == NULL ) {
    		perror("Errore in apertura del file");
  	}
  	
	for(int j = 0; j < N; j++){
		for(int i = 0; i <= M; i++){
			fprintf(fd, "%.20f \n", X[j][i]);
		}
	}
	
  	fclose(fd);
	return 0;
}
