function A = b5d(n)
    I = ones(n, 1);
    D = spdiags(-I, 0, n, n);
    T = spdiags([I, -2*I, I], [-1, 0, 1], n, n);
    A = kron(D, T) + kron(T, D);
endfunction
