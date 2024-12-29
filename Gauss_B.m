function B = Gauss_B(A)
D = diag(diag(A)); U = triu(A, 1); L = tril(A, -1);
B = -(D+L)\U;
end