A = [1 0 0 1
    -1 1 0 1
    -1 -1 1 1
    -1 -1 -1 1];
X = [0.2; 0.3; 0.1; 0.9]; %столбец решений
B = A*X; %вычисление правой части
ans1 = A\B; %поиск решений решатором
N1 = norm(ans1 - X)/norm(X); %относительная норма1 ошибки

[Q,R] = qr(A); %qr-разложение
ans2 = R\(Q'*B); %решение через qr
N2 = norm(ans2 - X)/norm(X); %относительная норма2 ошибки
[N1, N2]

vN1 = [];
vN2 = [];
for size=5:50
    E_big = eye(size);
    A_big = [E_big(:,1:end-1) ones(size,1)];
    X_big = rand(size,1);
    for k=1:size
      for col = 1:size
        for row = 1:size
            if col < row
                A_big(row,col)=-1;
            end
        end
      end
    end

    B_big = A_big*X_big; %вычисление правой части
    ans1_big = A_big\B_big; %поиск решений решатором
    N1_big = norm(ans1_big - X_big)/norm(X_big); %относительная норма1 ошибки
    vN1 = [vN1; N1_big];

    [Q_big,R_big] = qr(A_big); %qr-разложение
    ans2_big = R_big\(Q_big'*B_big); %решение через qr
    N2_big = norm(ans2_big - X_big)/norm(X_big); %относительная норма2 ошибки
    vN2 = [vN2; N2_big];
end
vN1

figure;

semilogy(vN1)
hold on;

semilogy(vN2)
hold off;
