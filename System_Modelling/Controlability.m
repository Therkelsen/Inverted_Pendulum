%% Controllability
% Define plant
A = [0    1.0000         0         0;
         0  -19.7863   25.8810   -0.0544;
         0         0         0    1.0000;
         0   -9.2336    2.2578   -0.0047]


B = [

         0;
    3.9573;
         0;
    1.8467]


C = [B A*B A^2*B A^3*B]
Cinv = inv(C)
rank(C)

% Find T
s4 = [0 0 0 1]*(Cinv)
s3 = s4*A
s2 = s3*A
s1 = s2*A

Tinv = [s1;s2;s3;s4]
T = inv(Tinv)


% Controllable canonical form Ac and Bc
Ac = Tinv*A*T
Bc = Tinv*B


% åben sløjfe polynomie og poler
syms lambda
eq = det(lambda*eye(4)-A)
r = vpasolve(eq,lambda)

Fc = [19.71-50 -2.6671-875 194.3013-6250 0-15000]

% feedback gains ud fra polplacering
F = Fc*Tinv 

% state space med feedback
Acf = Ac+Bc*Fc

% Check poles of plant with controller
roots([1 50 875 6250 15000])