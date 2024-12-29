function [x,k,B] = jacobi_or_gauss(A,b,x0, epslion,flag)
U = triu(A, 1); L =tril(A,-1); D = diag(diag(A));
B_J = -D\(L+U); B_G = -(D+L)\U;
g_J = D\b; g_G = (D+L)\b;
k = 1;
switch flag
    case 1
        x1 = B_J*x0+g_J;
        while norm(x1-x0,'inf') > epslion 
            x0 = x1;
            x1 = B_J*x0+g_J;
            k = k+1;
        end
        B = B_J;
    case 2
        x1 = B_G*x0+g_G;
        while norm(x1-x0,'inf') > epslion
            x0 = x1;
            x1 = B_G*x0+g_G;
            k = k+1;
        end
        B = B_G;
end
x = x1;