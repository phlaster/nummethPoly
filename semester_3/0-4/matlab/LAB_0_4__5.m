% 5б
clc;
repeats = 10;
n_max = 35;
matrix_sizes = round(logspace(log10(10), log10(2000), n_max));
times = zeros(1, n_max);

for repeat = 1:repeats
for i = 1:n_max
    n = matrix_sizes(i);
    A = fix_cond(n, 10);
    b = randn(n, 1);
   
    tic;
    x = A \ b;
    times(i) = times(i) + toc;
end
end
times = times / repeats;

loglog(matrix_sizes, times, 'o');
xlabel('Размер матрицы');
ylabel('Время вычисления, сек.');
title("Время решения СЛАУ оператором '\\'");
legend("Усреднение по " + string(repeats))
grid on;
hold on;
