#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define TYPE double
#define dx 0.001


void inicializar(TYPE p[], TYPE rho[], TYPE v[], TYPE x[], int N);
void inicialUF(TYPE u1[],TYPE u2[],TYPE u3[],TYPE f1[],TYPE f2[],TYPE f3[], TYPE rho[], TYPE v[], TYPE p[], int N);
TYPE E(TYPE v,TYPE rho,TYPE p);
void LaxWendroff(TYPE u1[],TYPE u2[],TYPE u3[],TYPE f1[],TYPE f2[],TYPE f3[], int N);
TYPE P(TYPE u1,TYPE u2,TYPE u3);
void actualizaF(TYPE u1[],TYPE u2[],TYPE u3[],TYPE f1[],TYPE f2[],TYPE f3[], int N);

int main()
{
  int l = 1;
  int N = l/dx;

  TYPE p[N];
  TYPE rho[N];
  TYPE v[N];
  TYPE x[N];
  
  inicializar(p,rho,v,x,N);
  
  TYPE u1[N];
  TYPE u2[N];
  TYPE u3[N];
  TYPE f1[N];
  TYPE f2[N];
  TYPE f3[N];

  inicialUF(u1,u2,u3,f1,f2,f3,rho,v,p,N);

  int j;

  //for(j=0;j<1000;j++)
  while( u2[(int)(0.9/dx)]/u1[(int)(0.9/dx)] - u2[(int)(0.9/dx + 1)]/u1[(int)(0.9/dx + 1)] == 0 )
    {
      LaxWendroff(u1,u2,u3,f1,f2,f3,N);
      actualizaF(u1,u2,u3,f1,f2,f3,N);
    }
  int i;
  for(i=0;i<N;i++)
  {
     /* Imprime posicion, densidad, presion y velocidad */
    printf("%f %f %f %f\n",x[i],u1[i],P(u1[i],u2[i],u3[i]),u2[i]/u1[i]);
  }
  return 0;
  
}

/* Inicializa densidad, presion, velocidad y la  posicion  */
void inicializar(TYPE p[],TYPE rho[], TYPE v[],TYPE x[], int N)
{
  int Naux=0.5/dx;
  int i;
  for(i=0;i<N;i++)
    {
      if(i<=Naux)
	{
	  p[i]=1.0;
	  rho[i]=1.0;
	}
      else
	{
	  p[i]=0.1;
	  rho[i]=0.125;
	}
      x[i]=i*dx;
      v[i]=0;
    }
}

/* Inicializa U y F */
void inicialUF(TYPE u1[],TYPE u2[],TYPE u3[],TYPE f1[],TYPE f2[],TYPE f3[], TYPE rho[], TYPE v[], TYPE p[], int N)
{
  int i;
  for(i=0;i<N;i++)
    {
      u1[i]=rho[i];
      u2[i]=rho[i]*v[i];
      u3[i]=E(v[i],rho[i],p[i]);

      f1[i]=rho[i]*v[i];
      f2[i]=p[i]+rho[i]*v[i]*v[i];
      f3[i]=v[i]*(E(v[i],rho[i],p[i])+p[i]);
    } 
}

TYPE E(TYPE v,TYPE rho,TYPE p)
{
  TYPE g=0.4;
  return rho*p/(rho*g)+rho*(pow(v,2))/2;
}

TYPE P(TYPE u1,TYPE u2,TYPE u3)
{
  return 0.4*(u3-0.5*u2*u2/(u1));
}

/* Actualiza el vector F */

void actualizaF(TYPE u1[],TYPE u2[],TYPE u3[],TYPE f1[],TYPE f2[],TYPE f3[], int N)
{
  int i;
  TYPE p;
  for(i=0;i<N;i++)
    {
      p = P(u1[i],u2[i],u3[i]);
      f1[i]=u2[i];
      f2[i]=p+u2[i]*u2[i]/u1[i];
      f3[i]=(u2[i]/u1[i])*(u3[i]+p);
    }
}

void LaxWendroff(TYPE u1[],TYPE u2[],TYPE u3[],TYPE f1[],TYPE f2[],TYPE f3[], int N)
{
  TYPE TU1[N];
  TYPE TU2[N];
  TYPE TU3[N]; 

  int r;
  for(r=0;r<N;r++)
    {
       TU1[r]=u1[r];
       TU2[r]=u2[r];
       TU3[r]=u3[r];
    }
  
  TYPE U1H1;
  TYPE U2H1;
  TYPE U3H1;
  TYPE U1H2;
  TYPE U2H2;
  TYPE U3H2;

  TYPE F1H1;
  TYPE F2H1;
  TYPE F3H1;
  TYPE F1H2;
  TYPE F2H2;
  TYPE F3H2;

  TYPE dt=0.0001;
  int i;
  for (i=1;i<N-1;i++)
    {  
      U1H1=0.5*(TU1[i+1]+TU1[i])-0.5*(dt/dx)*(f1[i+1]-f1[i]);
      U2H1=0.5*(TU2[i+1]+TU2[i])-0.5*(dt/dx)*(f2[i+1]-f2[i]);
      U3H1=0.5*(TU3[i+1]+TU3[i])-0.5*(dt/dx)*(f3[i+1]-f3[i]);
      U1H2=0.5*(TU1[i]+TU1[i-1])-0.5*(dt/dx)*(f1[i]-f1[i-1]);
      U2H2=0.5*(TU2[i]+TU2[i-1])-0.5*(dt/dx)*(f2[i]-f2[i-1]);
      U3H2=0.5*(TU3[i]+TU3[i-1])-0.5*(dt/dx)*(f3[i]-f3[i-1]);

      F1H1=U2H1;
      F2H1=P(U1H1,U2H1,U3H1)+U2H1*U2H1/U1H1;
      F3H1=(U2H1/U1H1)*(U3H1+P(U1H1,U2H1,U3H1));
      F1H2=U2H2;
      F2H2=P(U1H2,U2H2,U3H2)+U2H2*U2H2/U1H2;
      F3H2=(U2H2/U1H2)*(U3H2+P(U1H2,U2H2,U3H2));

      u1[i]=TU1[i]-(dt/dx)*(F1H1-F1H2);
      u2[i]=TU2[i]-(dt/dx)*(F2H1-F2H2);
      u3[i]=TU3[i]-(dt/dx)*(F3H1-F3H2);
    }
}

//TYPE Umax()
//{
//  TYPE umax=0;
//  int i;
//  for (i=1;i<N-1;i++)
//    {
//      if(u1[i]>=u2[i] && u1[i]>=u3[i] && u1[i]>=umax)
//	{
//	  umax=u1[i];
//	  continue;
//	}
//    else if(u2[i]>=u1[i] && u2[i]>=u3[i] && u2[i]>=umax)
//	{
//	  umax=u2[i];
//	  continue;
//	}
//    else if()
//	{
//	  umax=u3[i];
//	  continue;
//	}
//  }
//return umax;
//}
