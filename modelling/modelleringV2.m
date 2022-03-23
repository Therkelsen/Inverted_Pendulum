syms t q1 q2 q1t q2t I1z L1 M1        % variables
L = (I1z*q1t^2)/2 + (L1^2*M1*q1t^2)/8
dLdqt= [diff(L,q1t), diff(L,q2t)]
dLdq = [diff(L,q1), diff(L,q2)]
syms q1_f(t) q2_f(t)  % functions
q1t_f(t)= diff(q1_f,t)
q2t_f(t)= diff(q2_f,t)
    % replace the variables with the functions
dLdq_f= subs(dLdq,{q1 q2 q1t q2t},{q1_f q2_f q1t_f q2t_f})
dLdqt_f= subs(dLdqt,{q1 q2 q1t q2t},{q1_f q2_f q1t_f q2t_f}) 
    % now we can solve the equation
dsolve(diff(dLdqt_f,t)-dLdq_f==0)