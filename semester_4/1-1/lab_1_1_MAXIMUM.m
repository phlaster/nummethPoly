clc;
% МАКСИМУМ (+1)
% 8. Перестроить зависимости заполнения и времени, создав усреднение
% результатов по нескольким матрицам одной размерности и одной
% плотности соответственно
iters = 100;
density = logspace(-2, 0, iters);
t1 = zeros(iters,1);
t2 = zeros(iters,1);
t3 = zeros(iters,1);

msize = 100;
N = 15;
pb = uiprogressdlg(...
    uifigure,'Title',"Усреднение по "+string(N)+" результатам");
for d=1:iters
    pb.Message = "Вычлисление для плотности: " + string(density(d));
    t_1 = zeros(N,1);
    t_2 = zeros(N,1);
    t_3 = zeros(N,1);
    
    for i=1:N
        M = sprand(msize, msize, density(d), 1) + 10*eye(msize);
        pb.Value = density(d);

        t_1(i) = timeit(@() chol(M));
        t_2(i) = timeit(@() srcm(M));
        t_3(i) = timeit(@() samd(M));
    end
    
    t1(d) = mean(t_1);
    t2(d) = mean(t_2);
    t3(d) = mean(t_3);
end
close(pb);
plot(density, t1, density, t2, density, t3)
grid on
legend('только Холецкий','symrcm','symamd','Location','south')
titletext = "Сравнение методов, усреднение по "+string(N)+" результатам";
title(titletext)
xlabel('Плотность заполнения матрицы')
ylabel('Время одного вычисления, с')