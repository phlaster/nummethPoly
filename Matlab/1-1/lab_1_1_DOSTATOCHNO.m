clc;
% ДОСТАТОЧНО (+1)
% 6. Построить зависимость заполнения множителей разложения Холецкого
% при применении chol к исходной матрице и к матрицам с переставленными
% строками и столбцами при помощи symrcm и symamd.
nonzero1 = [];
nonzero2 = [];
nonzero3 = [];
msize = 100:400;
for n=msize
    M = sprandsym(n, 0.06, 0.8, 1);
    
    Ch1 = chol(M);
    nonzero1 = [nonzero1 nnz(Ch1)];
    
    km = symrcm(M);
    nonzero2 = [nonzero2 nnz(chol(M(km, km)))];
    
    sm = symamd(M);
    nonzero3 = [nonzero3 nnz(chol(M(sm, sm)))];
end

plot(msize, nonzero1,msize, nonzero2,msize, nonzero3)
grid on
legend('Исходная','symrcm','symamd','Location','southeast')
xlabel ('Размерность матрицы')
ylabel('Отношение ненулевых элементов ко всем')


% 7. Что быстрее?
% Определить время работы алгоритма при решения СЛАУ с разложением
% Холецкого без перестановок и с перестановками. Создать цикл по
% плотности заполнения матрицы. Построить график времени от плотности
t1 = [];
t2 = [];
t3 = [];
% density = logspace(-2.5, 0, 150);
density = logspace(-5, 0, 150);
msize = 160;
pb = uiprogressdlg(uifigure,'Title','Подождите...')
for d=density
    M = sprandsym(msize, d, 0.8, 1);
    pb.Value = d;
    pb.Message = "Вычлисление для плотности: " + string(d);
    t1 = [t1 timeit(@() chol(M))];
    t2 = [t2 timeit(@() srcm(M))];
    t3 = [t3 timeit(@() samd(M))];
end
close(pb);
loglog(density, t1, density, t2, density, t3)
grid on
legend('только Холецкий','symrcm','symamd','Location','south')
title("Сравнение методов")
xlabel('Плотность заполнения матрицы')
ylabel('Время одного вычисления, с')