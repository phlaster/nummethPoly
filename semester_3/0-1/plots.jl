using CairoMakie
using LaTeXStrings
using CSV
using ColorSchemes

begin
	f1(x)   = x^4-x^3-2x^2+3x-3
	d1f1(x) = 4x^3-3x^2-4x+3
	d2f1(x) = 12x^2-6x-4

	f2(x)   = 3ℯ^x-5x-3
	d1f2(x) = 3ℯ^x-5
	d2f2(x) = 3ℯ^x
end

pic2 = let
	f = Figure(resolution = (1100, 500))
	ax1 = Axis(f[1:3, 1:3], title="Отрицательный интервал (1)")
	hlines!(0, color=:black)
	
	a1, b1 = -2.732, -1.693
	x1 = range(a1, b1, length=100)
	lines!(ax1, x1, f1.(x1), label = L"f_i", linewidth=6, color=:crimson)
	lines!(ax1, x1, d1f1.(x1), label = L"f_i'", linewidth=4, color=:tomato)
	lines!(ax1, x1, d2f1.(x1), label = L"f_i''", linewidth=2, color=:lightsalmon)


	a2, b2 = 0.5, 4
	ax2 = Axis(f[1:3, 4:6], title="Положительный интервал (1)")
	hlines!(0, color=:black)
	
	x2 = range(a2, b2, length=100)
	lines!(ax2, x2, f1.(x2), label = L"f_i", linewidth=6, color=:crimson)
	lines!(ax2, x2, d1f1.(x2), label = L"f_i'", linewidth=4, color=:tomato)
	lines!(ax2, x2, d2f1.(x2), label = L"f_i''", linewidth=2, color=:lightsalmon)

	a3, b3 = 0.6, 1.5
	ax3 = Axis(f[1:3, 7:9], title="Интервал для (2)")
	hlines!(0, color=:black)
	
	x3 = range(a3, b3, length=100)
	lines!(ax3, x3, f2.(x3), label = L"f_i", linewidth=6, color=:dodgerblue1)
	lines!(ax3, x3, d1f2.(x3), label = L"f_i'", linewidth=4, color=:skyblue2)
	lines!(ax3, x3, d2f2.(x3), label = L"f_i''", linewidth=2, color=:lightblue2)

	Legend(f[4, 3:4], ax1, "")
	Legend(f[4, 8], ax3, "")

	current_figure()
end
save("./tex/pic2.png", pic2)


let
	a1, b1 = 0.5, 3

	f = Figure(resolution = (800, 800))
	ax = Axis(f[1, 1], title=L"x^4-x^3-2x^2+3x-3=0")
	hlines!(0, color=:black)
	
	x1 = range(a1, b1, length=100)
	lines!(ax, x1, f1.(x1), label = L"f_1")
	lines!(ax, x1, d1f1.(x1), label = L"f_1'")
	lines!(ax, x1, d2f1.(x1), label = L"f_1''")


	axislegend(position = :lt, L">0")
	current_figure()
end


function draw_errors(fnames::Vector{String}; leg=fnames, title="", cols=ColorSchemes.tab10.colors)
    fig = Figure(resolution = (1400, 800))
    ax1 = Axis(fig[1,1],
        title="Количество шагов от погрешности, метод $title",
        xticks = [-16:2:0;],
		yticks = [0:2:100;],
        xtickformat = values -> [L"10^{%$(round(Int, value))}" for value in values],
		xlabel=L"\varepsilon",
		ylabel="n steps",
    )
    xlims!(-16, 1)
	ylims = ()
    for (i, name) in enumerate(fnames)
        f = CSV.File(name)
		ylims = extrema([ylims..., f.steps...])
        lines!(log10.(f.eps), f.steps,
			label = leg[i],
			color=cols[i],
			linewidth=5,
			alpha=0.6,
		)
    end
	ylims!(0, +(ylims...))
	axislegend(position = :lb)

	ax2 = Axis(fig[1, 2],
        title="Фактическая точность от погрешности, метод $title",
        xticks = [-16:2:0;],
		yticks = [-16:2:0;],
        xtickformat = values -> [L"10^{%$(round(Int, value))}" for value in values],
		ytickformat = values -> [L"10^{%$(round(Int, value))}" for value in values],
		xlabel=L"\varepsilon",
		ylabel="abs. err",
    )
    xlims!(-16, 1)

    for (i, name) in enumerate(fnames)
        f = CSV.File(name)
        lines!(log10.(f.eps), log10.(f.err),
			label = leg[i],
			color=cols[i],
			linewidth=5,
			alpha=0.6,
		)
    end
	axislegend(position = :lt)
    fig
end

mpd = draw_errors(
	["bisect1.csv", "bisect2.csv"],
	leg=[L"f_1(x)=x^4-x^3-2x^2+3x-3", L"f_2(x)=3e^x-5x-3"],
	title="биекции",
	cols=[:red, :royalblue]
)
save("./tex/mpd.png", mpd)

newton = draw_errors(
	["newton1.csv", "newton2.csv"],
	leg=[L"f_1(x)=x^4-x^3-2x^2+3x-3", L"f_2(x)=3e^x-5x-3"],
	title="Ньютона",
	cols=[:red, :royalblue]
)
save("./tex/newton.png", newt)


# draw_errors(
# 	["broken_newton(limits).csv", "broken_newton(limits+x0).csv"],
# 	leg=[L"f_1(x)~broken~limits", L"f_1(x)~broken~limits + x_0"],
# 	title="Ньютона"
# )
broken_newton = draw_errors(
	["newton1.csv", "newton2.csv", "broken_newton(limits).csv", "broken_newton(limits+x0).csv"],
	leg=[L"f_1(x)", L"f_2(x)", L"f_1(x)~broken~limits", L"f_1(x)~broken~limits + x_0"],
	title="Ньютона с нарушениями",
	cols=[:red, :royalblue, :pink2, :magenta]
)
save("./tex/broken_newton.png", broken_newton)


let 
    f = Figure()

Axis(f[1, 1], ytickformat = values -> ["$(value)kg" for value in values])
Axis(f[1, 3], ytickformat = values -> [L"10^{%$(value^2)}" for value in values])
Axis(f[1, 4], ytickformat = values -> [rich("$value", superscript("XY", color = :red))
                                       for value in values])

f
    
end
let
	f = Figure(resolution = (1100, 500))
	ax1 = Axis(f[1:3, 1:3], title="Отрицательный интервал (1)")
	hlines!(0, color=:black)
	
	a1, b1 = -(1+√2), -(1+ 1/3^(1/3))
	x1 = range(a1, b1, length=100)
	lines!(ax1, x1, f1.(x1), label = L"f_i", linewidth=6, color=:crimson)
	lines!(ax1, x1, d1f1.(x1), label = L"f_i'", linewidth=4, color=:tomato)
	lines!(ax1, x1, d2f1.(x1), label = L"f_i''", linewidth=2, color=:lightsalmon)


	a2, b2 = 0.5, 3
	ax2 = Axis(f[1:3, 4:6], title="Положительный интервал (1)")
	hlines!(0, color=:black)
	
	x2 = range(a2, b2, length=100)
	lines!(ax2, x2, f1.(x2), label = L"f_i", linewidth=6, color=:crimson)
	lines!(ax2, x2, d1f1.(x2), label = L"f_i'", linewidth=4, color=:tomato)
	lines!(ax2, x2, d2f1.(x2), label = L"f_i''", linewidth=2, color=:lightsalmon)

	a3, b3 = 0.6, 1.5
	ax3 = Axis(f[1:3, 7:9], title="Интервал для (2)")
	hlines!(0, color=:black)
	
	x3 = range(a3, b3, length=100)
	lines!(ax3, x3, f2.(x3), label = L"f_i", linewidth=6, color=:dodgerblue1)
	lines!(ax3, x3, d1f2.(x3), label = L"f_i'", linewidth=4, color=:skyblue2)
	lines!(ax3, x3, d2f2.(x3), label = L"f_i''", linewidth=2, color=:lightblue2)

	Legend(f[4, 3:4], ax1, "")
	Legend(f[4, 8], ax3, "")

	current_figure()
end