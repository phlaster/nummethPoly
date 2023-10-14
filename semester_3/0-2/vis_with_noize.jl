using Plots, LaTeXStrings
using Plots: plot, plot!
using CSV
using DataFrames

function plotNoize(name, S)
    C = name[2]=='1' ? :royalblue : :tomato
    df = CSV.read("CSVs/task5_"*name*"_.csv", DataFrame)
    
    plot(legendtitle=latexstring("$(count(!ismissing, df.x_uniformNodes_Noize))\\ узлов,\\ \\mathbf{K}=$(name[end-7:end-4])"); S...)

    plot!(
        df.x_uniformTab, df.y_uniformTab,
        label=S[:fname],
        color=C,
        line=:dash,
        ;S...
    )

    scatter!(
        df.x_uniformNodes_Noize, df.y_uniformNodes_Noize,
        label=latexstring("Равномерные\\ узлы\\ с\\ шумом"),
        color=:black,
        markersize=5,
        ;S...
    )

    scatter!(
        df.x_chebNodes_Noize, df.y_chebNodes_Noize,
        label=latexstring("Узлы\\ Чебышёва\\ с\\ шумом"),
        color=:green4,
        markersize=9,
        marker=:xcross,
    )
end

S_common = Dict(
    :size => (1000, 680),
    :background_color => :white, :legendfontsize => 13, :markerstrokewidth => 0,
    :grid => true,
    :gridstyle => :dash,
    :gridlinewidth => 2,
    :gridalpha => 0.2,
)
S1_L = merge(S_common, Dict(
    :xlims => (0.4, 2.9),
    :xticks => (0.5:0.25:3),
    :ylims => (1, 6),
    :legend => :topleft,
    :fname => L"f_1(x) = ctg(x) + x^2"
))
S2_L = merge(S_common, Dict(
    :xlims => (-2.5, 2.3),
    :xticks => (-2.5:0.25:2),
    :ylims => (-10, 25),
    :legend => :topright,
    :fname => L"f_2(x) = x^5 - 3.2x^3 + 2.5x^2 - 7x + 1.5"
))
S_err = merge(S_common, Dict(
    :yaxis => (L"\varepsilon", (1e-2, 1e10), :log),
    :xaxis => (L"N", (2, 37)),
    :legend => :topleft,
    :xticks => ((5:5:80)),
    :yticks => ((10.0 .^ (-16:1:10))),
    :linewidth => 2.5,
    # :alpha=>0.8,
))
apath = "~/Documents/Edu/4сем/Nummethods/Labs/3сем/0-2/tex/pics/"


for n = [(4, 27), (9, 28), (30, 29), (35, 30)]
    p = plotNoize("f1_$(n[1])_0.200000", S1_L)
    savefig(p, apath * "pic$(n[2]).png")
end

for n = [(4, 31), (6, 32), (15, 33), (30, 34)]
    p = plotNoize("f2_$(n[1])_0.200000", S2_L)
    savefig(p, apath * "pic$(n[2]).png")
end

function plot_with_noize(fname, S)
    df = CSV.read("CSVs/task5_$(fname)_.csv", DataFrame)
    plot(; S...)
    plot!(
        df.x_uniformTab, df.y_uniformTab,
        label=S[:fname],
        color=:black,
        line=:dash,
        linewidth=1.5,
        markersize=3,)

    scatter!(
        df.x_uniformNodes_Noize, df.y_uniformNodes_Noize,
        label=latexstring("Узлы\\ равномерной\\ решётки:\\ $(count(!ismissing, df.x_uniformNodes_Noize))\\ шт."),
        color=:black,
        # markeralpha=0.8,
        markersize=4,)

    plot!(
        df.x_uniformTab, df.y_uniformLagr_Noize,
        label=latexstring("\$L("*S[:fname][2:4]*")\$"),
        color=:crimson,
        # alpha=0.45,
        # line=:dash,
        linewidth=2.5,
        ; S...)

    plot!(
        df.x_uniformTab, df.y_uniformHerm_Noize,
        label=latexstring("\$H("*S[:fname][2:4]*")\$"),
        color=:royalblue,
        # alpha=0.45,
        # line=:dash,
        linewidth=2.5,
        ; S...)


        plot!(twinx(),
        df.x_uniformTab, df.err_uniformLagr_Noize,
        yaxis=(L"\varepsilon", (1e-16, 1e16), :log),
        xlims=S[:xlims],
        color=:crimson,
        alpha=0.45,
        # line=:dash,
        linewidth=1.2,
        label=false,
        yticks=((10.0 .^ (-4:2:8))),
        ; S...)
    plot!([0], [0],
        color=:crimson,
        alpha=0.45,
        # line=:dash,
        label=L"Профиль\ ошибки\ полинома\ Лагранжа"
    )


    plot!(twinx(),
        df.x_uniformTab, df.err_uniformHerm_Noize,
        yaxis=(L"\varepsilon", (1e-16, 1e16), :log),
        xlims=S[:xlims],
        color=:royalblue,
        alpha=0.45,
        # line=:dashdot,
        linewidth=1.2,
        label=false,
        yticks=((10.0 .^ (-4:2:8))),
        # ; S...
    )
    plot!([0], [0],
        color=:royalblue,
        alpha=0.45,
        # line=:dashdot,
        label=L"Профиль\ ошибки\ сплайна\ Эрмита"
    )
end

for n = [(4, 27), (9, 28), (20, 29), (35, 30)]
    p = plot_with_noize("f1_$(n[1])_0.200000", S1_L)
    savefig(p, apath * "pic$(n[2]).png")
end

for n = [(4, 31), (6, 32), (15, 33), (30, 34)]
    p = plot_with_noize("f2_$(n[1])_0.200000", S2_L)
    savefig(p, apath * "pic$(n[2]).png")
end

plot_with_noize("f1_4_0.200000", S1_L)

function get_max_errs(taskf="task2-3_f1", col=:err_uniformLagr)
    CSV.read("CSVs/$(taskf).csv", DataFrame)[!, col] |> maximum
end

get_max_errs.("task5_f1_35_0.200000_", [:err_uniformLagr_Noize, :err_uniformHerm_Noize])
get_max_errs.("task5_f2_30_0.200000_", [:err_uniformLagr_Noize, :err_uniformHerm_Noize])






function get_max_errs(taskf, col; lims)
    maxerrs = Float64[]
    
    for i=lims
        e = CSV.read("CSVs/$(taskf)_$(i)_0.200000_.csv", DataFrame)[!, col]
        push!(maxerrs, e |> maximum)
    end
    
    maxerrs
end

function plotErr6(fname_1, fname_2, lims, S)
    er1_Lagr = get_max_errs(fname_1, :err_uniformLagr_Noize, lims=lims)
    er2_Lagr = get_max_errs(fname_2, :err_uniformLagr_Noize, lims=lims)

    er1_Herm = get_max_errs(fname_1, :err_uniformHerm_Noize, lims=lims)
    er2_Herm = get_max_errs(fname_2, :err_uniformHerm_Noize, lims=lims)
    p = plot()

    plot!(
        lims, er1_Lagr,
        label=L"Трансцендентная\ -\ Лагранж",
        color=:royalblue1,
        ; S...
    )
    plot!(
        lims, er2_Lagr,
        label=L"Полином\ -\ Лагранж",
        color=:crimson,
        ; S...
    )

    plot!(
        lims, er1_Herm,
        label=L"Трансцендентная\ -\ Эрмит",
        color=:royalblue1,
        line=:dash
        ; S...
    )
    plot!(
        lims, er2_Herm,
        label=L"Полином\ -\ Эрмит",
        color=:crimson,
        line=:dash
        ; S...
    )

end

p35 = plotErr6("task5_f1", "task5_f2", 4:35, S_err)
savefig(p35, apath * "pic35.png")