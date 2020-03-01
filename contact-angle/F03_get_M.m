function [w2p] = F03_get_M(u,n,FACTORIAL_MAT)
u_plus = max(u,0);
sum = 0;
% Denominator = factorial(n-1);
Denominator = 6; % 3!
exponent = n-1;

for k=0:n
    A5 = max(0,u_plus-k);
    A6 = FACTORIAL_MAT(k+1,n);
    A7 = (-1) ^ k ;
    sum = sum + A7* A6 * (A5^exponent);
end
w2p = sum/Denominator;
end
