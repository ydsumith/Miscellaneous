function [FACTORIAL_MAT] = F01_factorial(n)
MAT = zeros(n+1,n);
for k=0:n
    for i=1:n
        if k <= i
            MAT(k+1,i)=  factorial(i)/(factorial(k)*factorial(i-k));
        end
    end
end
  FACTORIAL_MAT = MAT;
end