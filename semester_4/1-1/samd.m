function a = samd(M)
    km = symamd(M);
    a = chol(M(km, km));
end
