clc;
F = b5d(4);
spy(F)
full(F)

C = chol(F, "lower");
spy(C)
full(C)

iC = ichol(F);
spy(iC)
full(iC)

options = struct('type', 'ict', 'droptol', 0.3);
iC_t = ichol(F, options);
spy(iC_t)
full(iC_t)
