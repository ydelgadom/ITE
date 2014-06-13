#include<iostream>
#include<fstream>
#include<cmath>
#include<complex>

using namespace std;

double m=0.2, Pi=4*atan(1);
int Ns=4, Nt=8, V=Ns*Ns*Ns*Nt;

// rho and rhobar
double mu, r1, r2;

// auxiliary variables
double A, B, R;
  
int main()
{
  
  ofstream fout("obs_free.dat");

  for(int nmu=0; nmu<201; nmu++){

    complex<double> exN, itN, exChi, itChi;
  
    mu=nmu*0.01;
    r1=exp(mu)-1;
    r2=exp(-mu)-1;
    
    for(int k1=-Ns/2; k1<Ns/2; k1++){
    for(int k2=-Ns/2; k2<Ns/2; k2++){ 
    for(int k3=-Ns/2; k3<Ns/2; k3++){
    for(int k4=-Nt/2; k4<Nt/2; k4++){
  
      // momentum components
      double p1=2*Pi*k1/Ns, p2=2*Pi*k2/Ns, p3=2*Pi*k3/Ns, p4=2*Pi*(k4+0.5)/Nt; 
      double c=cos(p4), s=sin(p4), sm=sinh(mu), cm=cosh(mu);
  
      A = m+4-cos(p1)-cos(p2)-cos(p3);
      B = sin(p1)*sin(p1)+sin(p2)*sin(p2)+sin(p3)*sin(p3);
      R = A*A-2*A*c+B+1;
      
      // additional auxiliary variables
      complex<double> a1 ((A*c-1.)/R, A*s/R);
      complex<double> a2 ((A*c-1.)/R,-A*s/R);
      complex<double> chi1 (-2*A*c*cm, -2*A*s*sm);
      complex<double> chi2 (A*A-2*A*c*cm+B+1, -2*A*s*sm);
      complex<double> chi3 (-2*A*c*sm, -2*A*s*cm);
      complex<double> sumN, sumChi;
  
      // exact particle number density
      exN += -2*A*sm*((A*A+B+1)*c-2*A*cm)/
             (pow(A*A-2*A*c*cm+B+1,2)+pow(2*A*s*sm,2));
      
      // exact particle number susceptibility
      exChi += chi1/chi2-chi3*chi3/(chi2*chi2);
  
      // 10th order ITE particle number density
      for(int i=1; i<2; i++){
        sumN += pow(a1*r1+a2*r2,i-1)*((a1)*exp(mu)-a2*exp(-mu));
      }
      itN += sumN;
      
      // 10th order ITE particle number susceptibility
      for(int j=1; j<3; j++){
        sumChi += double(j-1)*pow(a1*r1+a2*r2,j-2)*pow(a1*exp(mu)-a2*exp(-mu),2)+
                  pow(a1*r1+a2*r2,j-1)*(a1*exp(mu)+a2*exp(-mu));
      }
      itChi += sumChi;
      
    }
    }
    }   
    }
  
    // outfile
    fout << mu << " "
         <<  2.*real(exN)/double(V) << " "
         << -2.*real(itN)/double(V) << " "
         <<  2.*real(exChi)/double(V) << " "
         << -2.*real(itChi)/double(V) << endl;

  }
  
  fout.close();

return 0;
}
