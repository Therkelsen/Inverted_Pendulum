clc
clear 
close all

b1 = 1
b2 = 0.5
r = 1
J = 1
m = 5
b3 = 1


A = [0  1   0   0;0 -b1/J  0  -b2/J; 0 0 0 1; 0 0 0 -b3/m]
B = [0 0; -r/J r/J; 0 0; 1/m 1/m]   
C = [1 0 0 0; 0 0 1 0]


Obs = [C; C*A; C*A*A; C*A*A*A]

Con = [B  A*B   A*A*B   A*A*A*B]

nul = [0, 0; 0, 0; 0, 0; 0, 0]
nul1 = [0, 0; 0, 0]


Ae = [A, nul;C, nul1]
Be = [B; nul1]
Ce = [C, nul1]

rank(Obs)

rank(Con)

p = [-3.5 -3.6 -50 -51 -52 -53];
K = place(Ae,-Be,p)

Fi = K(:, [5, 6])
F = K(:, 1:4)

nul3 = [0, 0, 0]

p1 = [-50, -40, -45, -47]
L = place(A', -C', p1)
L = L'

