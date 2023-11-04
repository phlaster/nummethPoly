using CairoMakie, CSV, ColorSchemes, DataFrames

function plot_baza_1()
    f1 = CSV.read("dot-fish-2/data_baza_0.100000.csv", DataFrame)
    f2 = CSV.read("dot-fish-2/data_baza_0.050000.csv", DataFrame)
    kwargs = (;
    xminorticksvisible = true,
    xminorgridvisible = true,
    yminorticksvisible = true,
    yminorgridvisible = true,)
    
    f = Figure(resolution = (1000, 600))
    ax1 = CairoMakie.Axis(f[1, 1];
        backgroundcolor = :gray90,
        xlabel = "x", ylabel = "y",
        # xticks=[1:0.5:3;],
        xminorticks = IntervalsBetween(5),
        # yminorticks = IntervalsBetween(2),
        kwargs...
	)
    xlims!(0.9, 3.1)
    lines!(
        f2[!, 1], f2[!, 3],
        color = :grey20,
        label = "точно",
        )
    scatter!(
        f1[!, 1], f1[!, 2],
        color = :royalblue,
        label = "h=0.1",
        markersize=15
    )
    scatter!(
        f2[!, 1], f2[!, 2],
        color = :crimson,
        label = "h=0.05",
        markersize=9
    )
    axislegend(position=:rc)

    ax2 = CairoMakie.Axis(f[2, 1];
        backgroundcolor = :gray90,
        xlabel = "x", ylabel = "eps",
        # xticks=[1:0.5:3;],
        # yticks=[0:0.005:0.01;],
        xminorticks = IntervalsBetween(5),
        yminorticks = IntervalsBetween(5),
        kwargs...
    )
    xlims!(0.9, 3.1)
    scatterlines!(
        f1.x, f1[!, 4],
        color = :royalblue,
        label = "h=0.1",
        markersize=8
        )
    scatterlines!(
        f2.x, f2[!, 4],
        color = :crimson,
        label = "h=0.05",
        markersize=8
    )
    ylims!(0, 0.012)
    axislegend(position=:rc)
    # Legend(f[1, 2], ax1)
    # Legend(f[2, 2], ax2)
    f
end
plot_baza_1()


function plot_baza_2()
	f = CSV.read("dot-fish-2/data_baza_8dots.csv", DataFrame)
	fig = Figure(resolution = (1000, 600))
	tickformat(values) = [L"10^{%$(round(Int, log10(value)))}" for value in values]

	ax1 = Axis(
		fig[1:3, 1],
		xlabel = "Шаг разбиения",
		ylabel = "максимальная ошибка на интервале",
		xscale=log10,
		yscale = log10,
        xticks = 10. .^ [-10:0;],
        yticks = 10. .^ [-16:0;],
		xtickformat = tickformat,
        ytickformat = tickformat,
        backgroundcolor = :gray90,
	)
	xlims!(1e-9,1e0)
	ylims!(1e-14,1e0)
    lines!(ax1, 
        f[!, 1], f[!, 2],
        color=:grey,
        label = L"h^2"
    )
    scatterlines!(ax1, 
        f[!, 1], f[!, 3],
        color = :crimson,
        label = L"\sup(\Delta y)"
    )

    axislegend(position=:lt)
    current_figure()
end
plot_baza_2()


function plot_min_dost()
    f1 = CSV.read("dot-fish-2/data_min.csv", DataFrame)
    # f2 = DataFrame(
    #     eps= 10. .^ [-14:-1;],
    #     N= rand(9.2:0.1:10.8, 14) .^ [7:-0.5:0.5;],
    #     max_err = rand(10:0.1:11, 14) .^ [-14:-1;],
    # )
    f2 = CSV.read("dot-fish-2/data_dost.csv", DataFrame)
   
    fig = Figure(resolution = (1000, 600))
	tickformat(values) = [L"10^{%$(round(Int, log10(value)))}" for value in values]

	ax1 = Axis(
		fig[1, 1],
		# xlabel = "Шаг разбиения",
		ylabel = "N",
		xscale=log10,
		yscale = log10,
        xticks = 10. .^ [-16:0;],
        yticks = 10. .^ [0:8;],
		xtickformat = tickformat,
        ytickformat = tickformat,
        backgroundcolor = :gray90,
	)
	xlims!(1e-16, 1e0)
	# ylims!(1e-8,1e0)
    scatterlines!(ax1, 
        f1[!, 1], f1[!, 2],
        color=:crimson,
        label = "удвоение шагов"
    )
    scatterlines!(ax1, 
        f2[!, 1], f2[!, 2],
        color=:royalblue,
        label = "адаптивный"
    )
    axislegend(position=:lc)

    ax2 = Axis(
		fig[2, 1],
		xlabel = "Точность ε",
		ylabel = "Максимальная ошибка",
		xscale=log10,
		yscale = log10,
        xticks = 10. .^ [-16:0;],
        yticks = 10. .^ [-16:2:0;],
		xtickformat = tickformat,
        ytickformat = tickformat,
        backgroundcolor = :gray90,
	)
    xlims!(1e-16, 1e0)
	# ylims!(1e-8,1e0)
    lines!(ax2, [1e-15,1e0],[1e-15,1e0],
        color = :black,
        alpha = 0.6,
        label = "eps=err"
    )
    scatterlines!(ax2, 
        f1[!, 1], f1[!, 3],
        color=:crimson,
        # label = "удвоение шагов"
    )
    scatterlines!(ax2, 
        f2[!, 1], f2[!, 3],
        color=:royalblue,
        # label = "адаптивный"
    )
    axislegend(position=:lt)
    current_figure()
end
plot_min_dost()





y(x) = (3x^2*(x-1))^(1/3)
f(x,y) = x/y^2 + y/x
η(xi, yi, h) = f( xi + h/2,    yi + h/2 * f(xi, yi) )

rnd(x) = round.(x, digits=5)

function iter(x0, y0, h)
    @show e = η(x0, y0, h) |> rnd
    @show Δy = h * e |> rnd
    @show x1, y1 = x0 + h, y0 + Δy |> rnd
end

x0 = 1.1
y0 = y(x0)
h0 = 1/2
iter(x0, y0, h0)

h1 = 1/4
iter(x0, y0, h1)
iter(1.35, 1.20082, h1)


h2 = 1/8
iter(x0, y0, h2)
iter(1.225, 0.99281, h2)
iter(1.35, 1.2297, h2)
iter(1.475, 1.44647, h2)

1/3 * abs(1.57445 - 1.6236)
1/3 * abs(1.6236 - 1.65244)



hIII=1/8
x1, y1, = iter(x0, y0, hIII)
iter(x1, y1, hIII)
1/3 * abs(1.20082 - 1.2297)



hIV=1/16
x1, y1, = iter(x0, y0, hIV)
iter(x1, y1, hIV)
1/3 * abs(1.00151 - 0.99281)
