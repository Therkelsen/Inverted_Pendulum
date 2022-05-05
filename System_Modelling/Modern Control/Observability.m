%% Observability
clear
A = [   0   1   0   0
        0  -1.152  0.2817       0
        0       0       0       1
        0   2.468  -3.229       0   ]

B = [   0
        0.2304
        0
        -0.4937 ]
 
  C = [ 1   0   0   0
        0   0   1   0   ]

  D = [ 0
        0   ]

% Observability matrix
O = [   C;
        C*A;
        C*A^2;
        C*A^3   ]

obsv(A,C)

% Is observeable if full rank (rank = 4)
rank(O)


L = place(A, B, [-5, -10, -15, -20])