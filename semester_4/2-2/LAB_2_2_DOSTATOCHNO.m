% ДОСТАТОЧНО (+1) QR-алгоритм с предварительным приведением матрицы к форме
% Хессенберга
% – Построенную матрицу привести к форме Хессенберга ( hess() )
% – Добавить на графики линии для приведенной матрицы
% – Построить зависимость времени от размера матрицы для варианта приведенной
% матрицы и не приведенной


% v = [1:N];
% D = diag(v) + tril(rand(N), -1);
% [Q, R] = qr(rand(N));
% 
% A = Q' * D * Q;

N = 100;
M_size = 1:N;
time_no_hess = zeros(1,N);
time_hess = zeros(1,N);
for i=M_size
%       Создаём
    v = 1:i+9;
    D = diag(v) + tril(rand(i+9), -1);
    [Q, R] = qr(rand(i+9));
    A = Q' * D * Q;
      A1 = A;
  
%     Без приведения
    tic;
    while max(abs(tril(A1, -1)), [], "all") > 1e-2
        [Q,R] = qr(A1);
        A1 = R*Q;
    end
    time_no_hess(i) = toc;

%     С приведением
    A2 = A;
    tic;
        H = hess(A2);
    while max(abs(tril(H, -1)), [], "all") > 1e-2
        [Q,R] = qr(H);
        H = R*Q;
    end
    time_hess(i) = toc;    
end

semilogy(M_size, time_no_hess)
grid on
hold on
semilogy( M_size, time_hess)
hold on
legend("no Hess", "Hess")
title("Сравнение")
xlabel('M size')
ylabel('Time')

