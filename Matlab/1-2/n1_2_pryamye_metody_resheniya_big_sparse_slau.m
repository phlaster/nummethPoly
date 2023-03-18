clc
rng default %Задайте точку отсчета случайных чисел

A = rand(10,10);
A = A*A'; %случайную полностью заполненную симметричную ПО матрицу
X_t = ones(10,1); %вектор решений
C = A*X_t; % вычислите правую часть

X0 = pcg(A,C); %метод не сошелся за 10 итераций с точностью 1е-6
er = norm(C-A*X0)/norm(C) %проверка точности по формуле

X1 = pcg(A,C,1e-7,15);
er1 = norm(C-A*X1)/norm(C) %проверка точности по формуле


clc
tol = 1e-7; maxit = 15;
[x2,flag,rr2,it2,rv2] = pcg(A,C,tol,maxit)
% flag – описание сходимости метода
%    0 – метод сошелся с заданными входными параметрами
%    1 – число итераций = maxit, точность не достигнута
% rr2 – текущая достигнутая точность по формуле
% it2 – номер итерации для вычисленного x
% rv2 – вектор норм невязок на каждой итерации

clc
[x3] = pcg(A,C,1e-6,[],[],[],X0)
pogr(A, C, x3)


clc; clear;
% MINIMUM 1
B = sprandsym(10,0.5, 0.9, 1); spy(B);
X_t1 = ones(10,1);
C1 = B*X_t1;

X01 = pcg(B,C1);
err0 = pogr(B, C1, X01)
 
X11 = pcg(B,C1,1e-7,20);
err1 = pogr(B, C1, X11)



clc; clear;
% MINIMUM 2
R1 = rand(10,10);
X_t2 = ones(10,1);
C2 = R1*X_t2;

tol = 1e-7; maxit = 15;
[x2,flag,rr2,it2,rv2] = pcg(R1,C2,tol,maxit) % flag 4

% % X_t22 = ones(100,1);
% % R2 = sprandsym(100,0.2, 0.001, 1)
% % C22 = R2*X_t22;
% % [x2_,flag_,rr2_,it2_,rv2_] = pcg(R2,C22,tol,maxit) % flag 4


clc;
% DOST
N = 1000;
[Q, R] = qr(rand(N));
D = diag(1:N);
A = Q'*D*Q;
xt = ones(N,1);
Ca = A*xt;


B = sprandsym(1000,0.3,1e-3,1);
Cb = B*xt;

[xA,flagA,relresA,iterA] = pcg(A, Ca, 1e-7)
[xB,flagB,relresB,iterB] = pcg(B, Cb, 1e-7)

[xA200,flagA200,relresA200,iterA200] = pcg(A, Ca, 1e-7, 200)
[xB200,flagB200,relresB200,iterB200] = pcg(B, Cb, 1e-7, 200)

Ds = []
logspace(0,10,100)
%ФИГНЯЯ
for ObCh=logspace(0,10)
    Dn = diag(linspace(1, ObCh, 1000));
    xd = ones(1000,1);
    Cd = Dn*xd;
    [xD,flagD,relresD,iterD] = pcg(Dn, Cd, 1e-7, 300);
    Ds = [Ds, iterD]
%     B = sprandsym(1000,ObCh,1e-3,1);
%     [xB,flagB,relresB,iterB] = pcg(B, Cb, 1e-7);
%     Bs = [Bs, iterB];
end

Ds