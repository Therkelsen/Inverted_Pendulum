%% Controllability
% Define plant
A = [0    1.0000         0         0
         0  -19.7863   25.8810   -0.0544
         0         0         0    1.0000
         0   -9.2336    2.2578   -0.0047]


B = [    0
    3.9573
         0
    1.8467  ]

% Controlability
C = [B A*B A^2*B A^3*B]
rank(C)     % Controlable when full rank

% Find T
s4 = [0 0 0 1]*(inv(C))
s3 = s4*A
s2 = s3*A
s1 = s2*A

Tinv = [s1;s2;s3;s4]
T = inv(Tinv)


% Controllable canonical form Ac and Bc
Ac = Tinv*A*T
Bc = Tinv*B


% åben sløjfe polynomie
syms lambda
openLoopPoly = expand(simplify(det(lambda*eye(4)-A)))
openLoopPolyCoefs = [19.7910 2.6671  194.3013    0]

% lukket sløjfe polunomie
expand((lambda + 5)*(lambda + 10)*(lambda + 15)*(lambda + 20))
closedLoopPolyCoefs = [50  875 6250 15000]


% feedback gains ud fra polplacering
Fc = closedLoopPolyCoefs - openLoopPolyCoefs;
F = Fc*Tinv 

% state space med feedback
Acf = Ac+Bc*Fc

% Check poles of plant with controller
roots([1 50 875 6250 15000])

