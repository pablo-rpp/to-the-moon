#Llamada a Parámetros
param N;
#Iteradores
set D := {1 .. 2};
set NN := {0 .. N};
set NI := {1 .. N};
param pi := 4*atan(1); #Pi

param G;
param roE;
param roM;
param mE;
param mM;
param L;
param rE{D};
#param rE2;
param rM{D};
#param rM2;
param thetaE := 3 *pi/2; #Ángulo de salida
param thetaM := pi/2; #Ángulo de entrada

param TN;
param v0;
param vN;
param h:= TN/N;
#Velocidad inicial
param v01:= v0*cos(thetaE);
param v02:= v0*sin(thetaE);
#Velocidad final
param vN1:= vN*cos(pi - thetaM);
param vN2:= vN*sin(pi - thetaM);

#Posición inicial
param s01:= roE*cos(thetaE);
param s02:= roE*sin(thetaE);
#Posición final
param sN1:= L + roM*cos(thetaM);
param sN2:= roM*sin(thetaM);



#Variables
var s {D, NN}; #Posición
var v {D, NN}; #Velocidad
var a {D, NN}; #Aceleración gravitacional
var Th {D, NN}; #Aceleración artificial

#Función objetivo
minimize empuje: sum {d in D, n in NI} Th[d, n]^2;

subject to init_vel1 : v[1, 0] = v01;
subject to init_vel2 : v[2, 0] = v02;

subject to end_vel1 : v[1, N] = vN1;
subject to end_vel2 : v[2, N] = vN2;

subject to init_pos1 : s[1, 0] = s01;
subject to init_pos2 : s[2, 0] = s02;

subject to end_pos1 : s[1, N] = sN1;
subject to end_pos2 : s[2, N] = sN2;

subject to acc {d in D, n in NN}:
  a[d, n] = -G*mE*(s[d, n]-rE[d])/(sum{dd in D} (s[dd, n] - rE[dd])^2)^(3/2)
            -G*mM*(s[d, n]-rM[d])/(sum{dd in D} (s[dd, n] - rM[dd])^2)^(3/2)
            + Th[d, n]; #Auxiliar. Suma de aceleraciones en momento n

subject to vel_def {d in D, n in NI}:
  s[d,n] = (v[d, n-1] + v[d, n])*h/2 + s[d, n-1]; #ED posición

subject to acc_def {d in D, n in NI}:
  v[d, n] = v[d, n-1] + (a[d, n-1] + a[d, n])*h/2; #ED aceleración
