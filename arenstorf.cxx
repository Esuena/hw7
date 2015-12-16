#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

void GiveMe_Kvector(double* k, double w1, double w2, double w3, double w4){
  double mu=0.012277471;
  double r=sqrt((w1+mu)*(w1+mu)+w2*w2);
  double s=sqrt((w1-1+mu)*(w1-1+mu)+w2*w2);
  
 k[0]= w3;
 k[1]= w4;
 k[2]= w1+2*w4-(1-mu)*(w1+mu)/(r*r*r)-mu*(w1-1+mu)/(s*s*s);
 k[3]= w2-2*w3-(1-mu)*w2/(r*r*r)-mu*w2/(s*s*s);
}

void Dorman_Prince(double* y5, double* y, double& dt){
  double y1_temp, y2_temp, y3_temp, y4_temp;
  double k1[4], k2[4], k3[4], k4[4], k5[4], k6[4], k7[4];
  
  GiveMe_Kvector(k1, y[0], y[1], y[2], y[3]);
  
  y1_temp= y[0] + dt*k1[0]/5;
  y2_temp= y[1] + dt*k1[1]/5;
  y3_temp= y[2] + dt*k1[2]/5;
  y4_temp= y[3] + dt*k1[3]/5;
  
  GiveMe_Kvector(k2, y1_temp, y2_temp, y3_temp, y4_temp);
  
  y1_temp= y[0] + dt*(3*k1[0]/40 + 9*k2[0]/40);
  y2_temp= y[1] + dt*(3*k1[1]/40 + 9*k2[1]/40);
  y3_temp= y[2] + dt*(3*k1[2]/40 + 9*k2[2]/40);
  y4_temp= y[3] + dt*(3*k1[3]/40 + 9*k2[3]/40);
  
  GiveMe_Kvector(k3, y1_temp, y2_temp, y3_temp, y4_temp);
  
  y1_temp= y[0] + dt*(44*k1[0]/45 - 56*k2[0]/15 + 32*k3[0]/9);
  y2_temp= y[1] + dt*(44*k1[1]/45 - 56*k2[1]/15 + 32*k3[1]/9);
  y3_temp= y[2] + dt*(44*k1[2]/45 - 56*k2[2]/15 + 32*k3[2]/9);
  y4_temp= y[3] + dt*(44*k1[3]/45 - 56*k2[3]/15 + 32*k3[3]/9);
  
  GiveMe_Kvector(k4, y1_temp, y2_temp, y3_temp, y4_temp);
  
  y1_temp= y[0] + dt*(19372*k1[0]/6561 - 25360*k2[0]/2187 + 64448*k3[0]/6561 - 212*k4[0]/729);
  y2_temp= y[1] + dt*(19372*k1[1]/6561 - 25360*k2[1]/2187 + 64448*k3[1]/6561 - 212*k4[1]/729);
  y3_temp= y[2] + dt*(19372*k1[2]/6561 - 25360*k2[2]/2187 + 64448*k3[2]/6561 - 212*k4[2]/729);
  y4_temp= y[3] + dt*(19372*k1[3]/6561 - 25360*k2[3]/2187 + 64448*k3[3]/6561 - 212*k4[3]/729);
  
  GiveMe_Kvector(k5, y1_temp, y2_temp, y3_temp, y4_temp);
  
  y1_temp= y[0] + dt*(9017*k1[0]/3168 - 355*k2[0]/33 + 46732*k3[0]/5247 - 49*k4[0]/176 - 5103*k5[0]/18656);
  y2_temp= y[1] + dt*(9017*k1[1]/3168 - 355*k2[1]/33 + 46732*k3[1]/5247 - 49*k4[1]/176 - 5103*k5[1]/18656);
  y3_temp= y[2] + dt*(9017*k1[2]/3168 - 355*k2[2]/33 + 46732*k3[2]/5247 - 49*k4[2]/176 - 5103*k5[2]/18656);
  y4_temp= y[3] + dt*(9017*k1[3]/3168 - 355*k2[3]/33 + 46732*k3[3]/5247 - 49*k4[3]/176 - 5103*k5[3]/18656);
  
  GiveMe_Kvector(k6, y1_temp, y2_temp, y3_temp, y4_temp);
  
  y1_temp= y[0] + dt*(35*k1[0]/384 + 500*k3[0]/1113 + 125*k4[0]/192 - 2187*k5[0]/6784 + 11*k6[0]/84);
  y2_temp= y[1] + dt*(35*k1[1]/384 + 500*k3[1]/1113 + 125*k4[1]/192 - 2187*k5[1]/6784 + 11*k6[1]/84);
  y3_temp= y[2] + dt*(35*k1[2]/384 + 500*k3[2]/1113 + 125*k4[2]/192 - 2187*k5[2]/6784 + 11*k6[2]/84);
  y4_temp= y[3] + dt*(35*k1[3]/384 + 500*k3[3]/1113 + 125*k4[3]/192 - 2187*k5[3]/6784 + 11*k6[3]/84);
  
  GiveMe_Kvector(k7, y1_temp, y2_temp, y3_temp, y4_temp);
  
 
  //we continue the integration using the result y[i] from the fourth-order scheme
  
  for(int i=0; i<4; i++){
  //5th order
     y5[i] = y[i] + dt*(k1[i]*35/384 + k3[i]*500/1113 + k4[i]*125/192 - k5[i]*2187/6784 + k6[i]*11/84);
//   y5 = y + dt*(k1[1]*35/384 + k3[1]*500/1113 + k4[1]*125/192 - k5[1]*2187/6784 + k6[1]*11/84);
//   dx5 = dx + dt*(k1[2]*35/384 + k3[2]*500/1113 + k4[2]*125/192 - k5[2]*2187/6784 + k6[2]*11/84);
//   dy5 = dy + dt*(k1[3]*35/384 + k3[3]*500/1113 + k4[3]*125/192 - k5[3]*2187/6784 + k6[3]*11/84);
  
  //4th order
     y[i] += dt*(k1[i]*5179/57600 + k3[i]*7571/16695 + k4[i]*393/640 - k5[i]*92097/339200 + k6[i]*187/2100 + k7[i]/40);
//   y += dt*(k1[1]*5179/57600 + k3[1]*7571/16695 + k4[1]*393/640 - k5[1]*92097/339200 + k6[1]*187/2100 + k7[1]/40);
//   dx += dt*(k1[2]*5179/57600 + k3[2]*7571/16695 + k4[2]*393/640 - k5[2]*92097/339200 + k6[2]*187/2100 + k7[2]/40);
//   dy += dt*(k1[3]*5179/57600 + k3[3]*7571/16695 + k4[3]*393/640 - k5[3]*92097/339200 + k6[3]*187/2100 + k7[3]/40);
  }
}

void stepcontrol(double a, double b, double c, double d, double& dt, double& TOL){
  double V1;
  double V2;
  double Vmax=0;
  
  if(a>b)
    V1=a;
  else
    V1=b;
  
  if(c>d)
    V2=c;
  else
    V2=d;
  
  if(V1>V2)
    Vmax=V1;
  else
    Vmax=V2;
  
//calculate new stepsize
  dt *= 0.9 * pow(TOL/Vmax, 0.2); //0.2 wegen 1/(p+1), also 1/5
}



int main(){
 
  double y4[4], y5[4];//fuer 4th order
  
  //initial conditions
  y4[0] = 0.994;
  y4[1] = 0;
  y4[2] = 0;
  y4[3] = -2.00158510637908;
  
  double dt = 1e-3;
  double TOL = 1e-5;
  double T = 17.065216560157;
  
  double t=0; //set time
  
  ofstream myfile("data.txt");
    
    myfile << t << "\t" << dt << "\t" << y4[0] << "\t" << y4[1]  << endl;
    
  while(t<T){
    
    Dorman_Prince(y5, y4, dt);
    
    t += dt;
    
    myfile << t << "\t" << dt << "\t" << y4[0] << "\t" << y4[1]  << endl;
    
    stepcontrol(abs(y4[0]-y5[0]), abs(y4[1]-y5[1]), abs(y4[2]-y5[2]), abs(y4[3]-y5[3]), dt, TOL);
  }
  
  myfile.close();
  
 return 0;
}