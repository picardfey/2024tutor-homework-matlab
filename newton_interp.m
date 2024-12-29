function [poly,interdiff] = newton_interp(x, y)
n = length(x);
interdiff = zeros(n, n);
interdiff(:, 1) = y;
for j=2:n
    for i=j:n
        interdiff(i,j) = (interdiff(i, j-1)-interdiff(i-1, j-1))/(x(i)-x(i-j+1));
    end
end
coeff = diag(interdiff);
syms X
poly = coeff(1);
term = 1;
for i=2:n
    term = term*(X-x(i-1));
    poly = poly+coeff(i)*term;
end
