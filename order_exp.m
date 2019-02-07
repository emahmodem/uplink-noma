function [x] = order_exp(y)
% This function gets an unordered vector drawn from exponential distribution 
% and return the ordered vector
% y ~ EXP()
% x is the ordered x1<x2<  ....

n = size(y,2);
x = zeros(1,n);
for i = 1:n
    for j=1:i
        x(i) = x(i) + y(j) / (n -j +1);
    end
end
end