using CairoMakie
using LaTeXStrings
using CSV, DataFrames
using ColorSchemes

using LinearAlgebra, Statistics

function pplot(fname)
	f = CSV.read(fname, DataFrame)
	fig = Figure(resolution = (1500, 700))
	tickformat(values) = [L"10^{%$(round(Int, log10(value)))}" for value in values]
	labls = [L"\lambda_1", L"\lambda_n"]
	colors = [:royalblue, :crimson]

	ax1 = Axis(
		fig[1:3, 1],
		title="Шаги алгоритма от выбранной точности",
		xlabel = "Заданная точность ε",
		ylabel = "N шагов",
		xscale=log10,
		yscale = log10,
		yticks = Int.([1e1, 5e1, 1e2, 5e2, 1e3, 5e3, 1e4]),
		xtickformat = tickformat,
	)
	xlims!(1e-16,1e0)
	
	for i in 1:2
		scatterlines!(
			ax1,
			f.eps,
			f[!, 2i],
			label = labls[i],
			color=colors[i],
			linewidth=4,
		)
	end
	Legend(fig[1, 2], ax1)

	ax2 = Axis(
		fig[1:3, 3],
		title="Фактическая ошибка от выбранной точности",
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
	for i in 1:2
		scatterlines!(
			ax2,
			f.eps,
			replace(f[!, 2i+1], 0=>missing),
			linewidth=4,
			color=colors[i],
		)
	end

	current_figure()
end

p1 = pplot("CSVs/n=10_cond=4.csv")
p2 = pplot("CSVs/n=10_cond=8.csv")
p3 = pplot("CSVs/n=10_cond=16.csv")
p4 = pplot("CSVs/n=10_cond=32.csv")

p5 = pplot("CSVs/n=15_cond=4.csv")
p6 = pplot("CSVs/n=20_cond=4.csv")

pplot("CSVs/n=15_cond=100.csv")

# _path = "tex/pics/conv_"
# save.(
# 	_path .* ["rand", "cond=500", "cond=50", "cond=5"] .* ".png",
# 	[p1,p2,p3,p4]
# )


v0 = rand(4)
v = v0./norm(v0)
H = I - 2v*v'
D = diagm([4:-1:1;])
A = H' * D * H
xax = 0:0.001:5
conds = [cond(A-x*I) for x in xax]
let
	fig = Figure(resolution = (1000, 500))
	tickformat(values) = [L"10^{%$(round(Int, log10(value)))}" for value in values]

	ax = Axis(
		fig[1:3, 1],
		title="Влияние сдвига на число обусловленности",
		xlabel = "Скаляр сдвига μ",
		ylabel = "cond",
		yscale = log10,
		xticks = [0:0.5:5;],
		yminorgridvisible=true
		# yticks = Int.([1e1, 5e1, 1e2, 5e2, 1e3, 5e3, 1e4]),
		# xtickformat = tickformat,
	)
	ylims!(1e-0,1e4)
	lines!(
		ax, xax, conds, label="Ч.о. после сдвига", color=:red, linewidth=2)
	lines!(ax, [0,5], [4,4], label="Исходное ч.о.", linewidth=2)

	Legend(fig[1, 2], ax)

	
	current_figure()
end
