% 4б
clc;

for n = [5, 10, 15]
    A = hilb(n);
    b = A * randn(n, 1);
    x = A \ b;
    
    disp(['Размер: ' num2str(n)]);
    disp(['Число обусловленности: ' num2str(cond(A))]);
    disp(['Норма фактической ошибки: ' num2str(norm(A*x - b))]);
    disp(['Норма невязки: ' num2str(norm(A*x - b))]);
    disp('-----------------------------');
end

