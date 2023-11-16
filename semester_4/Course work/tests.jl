using LinearAlgebra, Statistics

L = rand(-10:10, 5,5) |> LowerTriangular |> Matrix
U = rand(-10:10, 5,5) |> UpperTriangular |> Matrix
x = rand(1:9, 5)
bL = L*x
bU = U*x
L\b

w0 = [4,5,1,6];
w = w0 ./ norm(w0);
H = I - 2w*w';
D = diagm([5,3,2,1]);
A = H'*D*H;

A[1,4] = A[4,1] = 0

C = cholesky(Hermitian(A))
C.L' * C.U'

function chol(A)
    n = size(A, 1)  # Получаем размер матрицы A
    L = zeros(n, n)  # Создаем нулевую матрицу L того же размера
    
    for j in 1:n
        print("j = $j")
        @show _su = j == 0 ? [0] : [L[j, k]^2 for k in 1:j-1]
        @show su = sum(_su)
        @show s = A[j, j] - su
        @show L[j, j] = sqrt(s)
        println("DIAG:")
        display(round.([A L], digits=2))
        for i in j+1:n
            print("i = $i")
            @show _su = j == 0 ? [0] : [L[i, k]*L[j, k] for k in 1:j-1]
            @show su = sum( _su )
            @show s = A[i, j] - su
            @show L[i, j] = s / L[j, j]
        end
        println("REST:")
        display(round.([A L], digits=2))
    end
    return L
end

A = rand(0:1, 5,5); A = A*A' + diagm(2 .*ones(Int, 5))
A .*= A
L = zeros(5,5)

C = chol(A)

round.([A L], digits = 1)

j += 1
L[j, j] = sqrt(A[j, j] - sum( [L[j, k]^2 for k in 1:j-1] ))

for i in j+1:5
    L[i, j] = (A[i, j] - sum( [L[i, k]*L[j, k] for k in 1:j-1] )) / L[j, j]
end

using CairoMakie, CSV, ColorSchemes, DataFrames


function plot_baza_1()
    tickformat(values) = [L"10^{%$(round(Int, log10(value)))}" for value in values]


    f1 = CSV.read("CSVs/100x100_cond=500.csv", DataFrame)
    kwargs = (;
    xminorticksvisible = true,
    xminorgridvisible = true,
    yminorticksvisible = true,
    yminorgridvisible = true,)
    
    f = Figure(resolution = (1000, 600))
    ax1 = CairoMakie.Axis(f[1, 1];
        backgroundcolor = :gray90,
        xlabel = "порог", ylabel = "шаги",
        xticks = 10. .^ [-10:4;],
        xtickformat = tickformat,
        xminorticks = IntervalsBetween(10),
        yminorticks = IntervalsBetween(10),
        xscale=log10,
        kwargs...
	)
    scatterlines!(
        f1[!, 1], f1[!, 2],
        color = :royalblue,
        label = "pcg"
    )
    scatterlines!(
        f1[!, 1], fill(27, 19),
        color = :crimson,
        label = "cg"
    )
    axislegend(position=:rb)
    f
end
plot_baza_1()