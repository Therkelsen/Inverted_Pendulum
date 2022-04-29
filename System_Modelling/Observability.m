%% Observability
clear
A = [0    1.0000         0         0;
         0  -19.7863   25.8810   -0.0544;
         0         0         0    1.0000;
         0   -9.2336    2.2578   -0.0047]


B = [   0;
        3.9573;
        0;
        1.8467  ]

% Multiple output (our system)
%C = [   1   0   0   0;
%        0   0   1   0   ]

% Single output (not our system)
C = [   1   0   0   0  ]

% Observability matrix
O = [   C;
        C*A;
        C*A^2;
        C*A^3   ]

% Is observeable if full rank (rank = 4)
rank(O)

% Find Observable canonical form
t4 = inv(O)*[0;0;0;1];
t3 = A*t4;
t2 = A*t3;
t1 = A*t2;
T = [t1 t2 t3 t4]