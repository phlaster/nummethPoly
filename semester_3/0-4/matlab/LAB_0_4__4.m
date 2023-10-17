% 4б
clc;

for n = [5, 10, 15]
    A = hilb(n);
    b = randn(n, 1);
    x = A \ b;

    cond_num = cond(A);
    actual_error = norm(A*x - b);
    residual = norm(A*x - b) / norm(b);
    
    disp(['Размер: ' num2str(n)]);
    disp(['Число обусловленности: ' num2str(cond_num)]);
    disp(['Норма фактической ошибки: ' num2str(actual_error)]);
    disp(['Норма невязки: ' num2str(residual)]);
    disp('-----------------------------');
end
