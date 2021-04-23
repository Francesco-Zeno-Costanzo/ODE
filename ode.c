#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define o0 9.0

double f(double v){
    double x_dot = v;
    return x_dot;
}

double g(double x){
    double v_dot = - o0*x;
    return v_dot;
}

double** Eulero(double x0, double v0, int N, double tf){
	double dt=tf/(N*1.0);
	double **X;
	X=(double **)malloc(3*sizeof(double));
	for(int i=0; i<3; i++){
    		X[i] = (double*)calloc(N+1, sizeof(double) );
	}
	X[0][0]=x0;
	X[1][0]=v0;
	for(int i=0; i<N; i++){
		X[0][i + 1] = X[0][i] + dt*f(X[1][i]);
		X[1][i + 1] = X[1][i] + dt*g(X[0][i]);
		X[2][i + 1] = X[2][i] + dt;
	}
	return X;
}

double** Eulero_implicito(double x0, double v0, int N, double tf){
	double dt=tf/(N*1.0);
	double **X;
	X=(double **)malloc(3*sizeof(double));
	for(int i=0; i<3; i++){
    		X[i] = (double*)calloc(N+1, sizeof(double) );
	}
	X[0][0]=x0;
	X[1][0]=v0;
	for(int i=0; i<N; i++){
		X[1][i + 1] = (X[1][i] -o0*dt*X[0][i])/(1.0+o0*dt*dt);
		X[0][i + 1] = X[0][i] + dt*X[1][i + 1];
		X[2][i + 1] = X[2][i] + dt;
	}
	return X;
}

double** Eulero_semi_implicito(double x0, double v0, int N, double tf){
	double dt=tf/(N*1.0);
	double **X;
	X=(double **)malloc(3*sizeof(double));
	for(int i=0; i<3; i++){
    		X[i] = (double*)calloc(N+1, sizeof(double) );
	}
	X[0][0]=x0;
	X[1][0]=v0;
	for(int i=0; i<N; i++){
		X[1][i + 1] = X[1][i] + dt*g(X[0][i]);
		X[0][i + 1] = X[0][i] + dt*f(X[1][i+1]);
		X[2][i + 1] = X[2][i] + dt;
	}
	return X;
}

double** velocity_verlet(double x0, double v0, int N, double tf){
	double dt=tf/(N*1.0);
	double **X;
	X=(double **)malloc(3*sizeof(double));
	for(int i=0; i<3; i++){
    		X[i] = (double*)calloc(N+1, sizeof(double) );
	}
	X[0][0]=x0;
	X[1][0]=v0;
	for(int i=0; i<N; i++){
    		X[0][i + 1] = X[0][i] + (X[1][i])*dt + g(X[0][i])*(dt*dt/2.0);
    		X[1][i + 1] = X[1][i] + (g(X[0][i])+g(X[0][i+1]))*dt/2.0;
		X[2][i + 1] = X[2][i] + dt;
	}
	return X;
}
	
double** RK4(double x0, double v0, int N, double tf){
	double dt=tf/(N*1.0);
	double **X;
	X=(double **)malloc(3*sizeof(double));
	for(int i=0; i<3; i++){
    		X[i] = (double*)calloc(N+1, sizeof(double) );
	}
	X[0][0]=x0;
	X[1][0]=v0;
	double xk1, xk2, xk3, xk4, vk1, vk2, vk3, vk4;
	for(int i=0; i<N; i++){
		xk1=xk2=xk3=xk4=vk1=vk2=vk3=vk4=0;
		xk1 = f(X[1][i]);
		vk1 = g(X[0][i]);
		xk2 = f(X[1][i] + vk1*dt/2.0);
		vk2 = g(X[0][i] + xk1*dt/2.0);
		xk3 = f(X[1][i] + vk2*dt/2.0);
		vk3 = g(X[0][i] + xk2*dt/2.0);
		xk4 = f(X[1][i] + vk3*dt);
		vk4 = g(X[0][i] + xk3*dt);
		X[0][i + 1] = X[0][i] + (dt/6.0)*(xk1 + 2.0*xk2 + 2.0*xk3 + xk4);
		X[1][i + 1] = X[1][i] + (dt/6.0)*(vk1 + 2.0*vk2 + 2.0*vk3 + vk4);
		X[2][i + 1] = X[2][i] + dt;
	}
	return X;
}
   
double** punto_medio_impliicto(double x0, double v0, int N, double tf){
	double dt=tf/(N*1.0);
	double **X;
	X=(double **)malloc(3*sizeof(double));
	for(int i=0; i<3; i++){
    		X[i] = (double*)calloc(N+1, sizeof(double) );
	}
	X[0][0]=x0;
	X[1][0]=v0;
	for(int i=0; i<N; i++){
    		X[1][i + 1] = (X[1][i] - o0*dt*X[0][i] - (o0*dt*dt*X[1][i])/4.0)/(1.0+(o0*dt*dt)/4.0);
		X[0][i + 1] = X[0][i] + dt*(X[1][i] + X[1][i + 1])/2.0;
		X[2][i + 1] = X[2][i] + dt;
	}
	return X;
}

double** predictor_corrector(double x0, double v0, int N, double tf, int n){
	double dt=tf/(N*1.0);
	double **X;
	double xp0, vp0, xp1, vp1, xc0, vc0, xc1, vc1;
	X=(double **)malloc(3*sizeof(double));
	for(int i=0; i<3; i++){
    		X[i] = (double*)calloc(N+1, sizeof(double) );
	}
	X[0][0]=x0;
	X[1][0]=v0;
	for(int i=0; i<N; i++){
		//xp0=vp0=xp1=vp1=xc0=vc0=xc1=vc1=0;
		xp0 = f(X[1][i]);
		vp0 = g(X[0][i]);
		xp1 = X[0][i] + dt*xp0;
		vp1 = X[1][i] + dt*vp0;
		X[2][i + 1] = X[2][i] + dt;
		xc0=f(vp1);
		vc0=g(xp1);
		for(int j=0; j<n; j++){
			xc1 = X[0][i] + 0.5*dt*(xp0 +xc0);
			vc1 = X[1][i] + 0.5*dt*(vp0 +vc0);
			xc0 = f(vc1);
			vc0 = g(xc1);
		}
		X[0][i + 1] = X[0][i] + 0.5*dt*(xp0 +xc0);
		X[1][i + 1] = X[1][i] + 0.5*dt*(vp0 +vc0);
	}
	return X;
}
	
int main(void){
	int x, N, n;
	double tf, x0, v0;
	double **X;
	X=(double **)malloc(3*sizeof(double));
	for(int i=0; i<3; i++){
    		X[i] = (double*)calloc(N+1, sizeof(double) );
	}
		
	printf("Per il metodo di eulero digitare 1 \n");
	printf("Per il metodo di eulero implicito digitare 2 \n");
	printf("Per il metodo di eulero semi implicito digitare 3 \n");
	printf("Per il metodo velocity verlet digitare 4 \n");
	printf("Per il metodo rk4 digitare 5 \n");
	printf("Per il metodo punto medio implicito digitare 6 \n");
	printf("Per il metodo predizione correzione digitare 7 \n");
	scanf("%d",&x);
	
	printf("ditanza:");
	scanf("%lf",&tf);
	printf("# di passi:");
	scanf("%d",&N);
	printf("x(t=0):");
	scanf("%lf",&x0);
	printf("v(t=0):");
	scanf("%lf",&v0);
	
	if (x==1) X=Eulero(x0, v0, N, tf);
	if (x==2) X=Eulero_implicito(x0, v0, N, tf);
	if (x==3) X=Eulero_semi_implicito(x0, v0, N, tf);
	if (x==4) X=velocity_verlet(x0, v0, N, tf);
	if (x==5) X=RK4(x0, v0, N, tf);
	if (x==6) X=punto_medio_impliicto(x0, v0, N, tf);
	if (x==7){
		printf("quante colte applicare la correzione? \n");
		scanf("%d",&n);
		X=predictor_corrector(x0, v0, N, tf, n);
	}
	
	FILE *fd;
	fd=fopen("scrivi.txt", "w");
	if( fd==NULL ) {
    		perror("Errore in apertura del file");
  	}
  	
	for(int j=0; j<=N; j++){
		fprintf(fd, "%.20f \t %.20f  \t %.20f \n", X[2][j], X[0][j], X[1][j]);
	}
	
  	fclose(fd);
	return 0;
}
