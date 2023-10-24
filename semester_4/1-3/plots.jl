using CairoMakie
using LaTeXStrings
using CSV, DataFrames
using ColorSchemes


function pplot(fname)
	f = CSV.read(fname, DataFrame)
	fig = Figure(resolution = (1500, 700))
	tickformat(values) = [L"10^{%$(round(Int, log10(value)))}" for value in values]

	ax1 = Axis(
		fig[1:3, 1],
		title="Шаги алгоритма от выбранной точности, усреднение по 10",
		xlabel = "Заданная точность ε",
		ylabel = "N шагов",
		xscale=log10,
		yscale = log10,
		yticks = [10, 30, 50, 100, 200, 500, 1000, 2000],
		xtickformat = tickformat,
	)
	xlims!(1e-16,1e0)
	
	for i in 1:6
		scatterlines!(
			ax1,
			f.eps,
			f[!, 2i],
			label = names(f)[2i][3:end],
			linewidth=4,
		)
	end
	Legend(fig[1, 2], ax1, "Размер")

	ax2 = Axis(
		fig[1:3, 3],
		title="Фактическая ошибка от выбранной точности (худшая из 10)",
		xlabel = "Заданная точность ε",
		ylabel = "Норма фактической ошибки",
		xscale=log10,
		yscale = log10,
		xticks = 10. .^[-16:2:0;],
		xtickformat = tickformat,
		yticks = 10. .^[-16:2:0;],
		ytickformat = tickformat,
	)
	xlims!(1e-16,1e0)
	lines!(
		[1e-16,1e0],[1e-16,1e0],color=:black, linestyle =:dash, )
	for i in 1:6
		scatterlines!(
			ax2,
			f.eps,
			f[!, 2i+1],
			linewidth=4,
		)
	end

	current_figure()
end

p1 = pplot("CSVs/pcg_convergence_randsym.csv")
p2 = pplot("CSVs/pcg_convergence_cond=500.csv")
p3 = pplot("CSVs/pcg_convergence_cond=50.csv")
p4 = pplot("CSVs/pcg_convergence_cond=5.csv")


_path = "tex/pics/conv_"
save.(
	_path .* ["rand", "cond=500", "cond=50", "cond=5"] .* ".png",
	[p1,p2,p3,p4]
)

let
	fig = Figure(resolution = (1000, 400))

	ax1 = Axis(
		fig[1, 1],
		title="Усреднение по 100",
		xscale=log10,
		yscale = log10,
		xticks = vcat(10, 20:20:100, 200:200:1000, 1500),
		yticks = 10. .^ [-15:2:0;],
		ytickformat = values -> [L"10^{%$(round(Int, log10(value)))}" for value in values],
	)
	
	f = CSV.read("CSVs/multi_x100.csv", DataFrame)
	for (i, n)=enumerate([10,30,50,100,200])
		scatterlines!(
			ax1,
			f[!, i+1],
			f[!, 1],
			label = "$n",
			linewidth=4,
		)
	end

	Legend(fig[1, 2], ax1, "Размер")

	current_figure()
end