function a = srcm(M)
    km = symrcm(M);
    a = chol(M(km, km));
end
