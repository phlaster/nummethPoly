clc;
% ДОСТАТОЧНО (+1)
% Исследовать зависимость количества итераций от числа обусловленности
% матрицы при одной и той же заданной точности
% 1. Создайте матрицы А и В ( плотность 0,3) размерностью 1000х1000
% с числом обусловленности 1000, единичный вектор решений xt и вычислите
% правую часть Са и Сb
N = 1000;
[Q, R] = qr(rand(N));
D = diag(1:N);
A = Q'*D*Q;
xt = ones(N,1);
Ca = A*xt;

B = sprandsym(1000,0.3,1e-3,1);
Cb = B*xt;

% 2. Примените pcg с точностью 1e-7 и выходными параметрами flag и iter.
[xA,flagA,relresA,iterA] = pcg(A, Ca, 1e-7);
[xB,flagB,relresB,iterB] = pcg(B, Cb, 1e-7);

% 3. Увеличьте число итераций до 200. Убедитесь, что метод сошелся для
% обеих матриц. Если нет – увеличьте максимальное число итераций
[xA200,flagA200,relresA200,iterA200] = pcg(A, Ca, 1e-7, 200);
[xB200,flagB200,relresB200,iterB200] = pcg(B, Cb, 1e-7, 200);

% 4. Создайте цикл для изменения числа обусловленности матриц от 1
% до 1е10 и постройте график зависимости числа итераций от числа
% обусловленности для обеих матриц
n = 200;
obus = logspace(0,10, n);
itersA = zeros(n,1);
itersB = zeros(n,1);
msize = 200;
maxiter = 4000;
pb = uiprogressdlg(uifigure,'Title','Подождите...');
for i=1:n
    pb.Value = i/n;
    xt = ones(msize,1);
    A = randcond(msize, obus(i));
    B = sprandsym(msize, 0.3, 1/obus(i), 1);
    Ca = A*xt;
    Cb = B*xt;
    
    [xA,flagA,relresA,iterA] = pcg(A, Ca, 1e-7, maxiter);
    [xB,flagB,relresB,iterB] = pcg(B, Cb, 1e-7, maxiter);
    
    itersA(i) = iterA;
    itersB(i) = iterB;
end
close(pb);
loglog(obus, itersA, obus, itersB);
grid on;
legend('Вручную (заполнение 100%)','sprandsym (заполнение 30%)',...
    'Location','northwest');
xlabel('Число обусловленности матрицы');
ylabel('Количество шагов для точности 1e-7');
