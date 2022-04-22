clc, clear, close all
%  
%                    3.957 s^6 + 78.24 s^5 + 26.69 s^4 + 1538 s^3 - 119.5 s^2 + 7551 s + 5.165e-12
%    1:  ------------------------------------------------------------------------------------------------------
%        s^8 + 39.58 s^7 + 386.4 s^6 + 283.1 s^5 + 7698 s^4 - 1036 s^3 + 3.775e04 s^2 - 3.065e-12 s - 9.241e-27
%  
%               1.847 s^4 + 36.55 s^3 - 4.924 s^2 + 358.8 s + 1.091e-11
%    2:  ----------------------------------------------------------------------
%        s^6 + 39.58 s^5 + 386.4 s^4 + 283.1 s^3 + 7698 s^2 - 1036 s + 3.775e04
% 

tfC = tf([3.957 78.24 26.69 1538 -119.5 7551  0], [1 39.58 386.4  283.1  7698  -1036 37750  0  0])
tfP = tf([1.847 36.55 -4.924 358.8 0], [1 39.58 386.4  283.1  7698  -1036   0])


[numC, denC] = tfdata(tfC)
[numP, denP] = tfdata(tfP)

s = tf('s')

t = 0:0.01:10;
u = zeros(size(t));
u(t == 1) = 1;

figure()
lsim(tfP, u, t)

K = 3.02 + 0.0933*(1/s) + 21.7*(s/(0.00165*s + 1))


