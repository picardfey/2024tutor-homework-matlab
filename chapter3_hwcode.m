%% 1
A = [1 -2 2
    -1 1 -1
    -2 -2 1];
b = [3 -2 -3]';
res = norm(A,2); % 1,inf
B = A*b;
disp(norm(B,'inf'))
%% simple iter 4 
clc,clear
% A = [5 -2 2
%     -1 5 -1
%     -2 -2 5];
% b = [5, 3, 1]';

% A = [2 1 1
%     1 1 1
%     1 1 2];
% b = [0 1 2]';

A = [1 -2 2
    -1 1 -1
    -2 -2 1];
b = [3 -2 -3]';

% A = [-1 5 -1
%     -2 -2 5
%     5 -2 2];
% b = [-3 5 9]';
% A = [5 -1 -1
%     -2 5 -2
%     -2 2 5];
x0 = zeros(3, 1);epslion = 5e-4; 
[x_jacobi1,k_jacobi1, B_J] = jacobi_or_gauss(A, b, x0, epslion, 1);
[x_gauss1, k_gauss1, B_G] = jacobi_or_gauss(A, b, x0, epslion, 2);
disp(['jacobi迭代解是',' gauss迭代解是'])
disp([x_jacobi1,x_gauss1])
disp(['jacbobi迭代次数是',num2str(k_jacobi1),'  gauss迭代次数是',num2str(k_gauss1)])
%% Richardson迭代 8 
clc,clear
A = [-8 6
    -4 2];
b = ones(2 ,1);
a = -1/3; k = 1;epslion = 5e-4;
x0 = [0.2 0.3]';x1 = x0 - a*(A*x0-b);
while norm(x1-x0,1)>epslion && k<30
    x0 = x1;
    x1 = x0 - a*(A*x0-b);
    k = k+1;
end
disp(x1)
%% 10 LU分解
clc,clear
A = [6 1 1
    1 6 1
    1 1 6];
b = [0 -15 -15]';
[L, U, D] = lu(A);
b = D\b;
%% cholesky分解
clc, clear
A = [16 4 8
    4 5 -4
    8 -4 22];
b = [1 1 1]';
[L, flag] = chol(A);
disp(A\b)