using Plots, ColorSchemes, LaTeXStrings
using CSV
using DataFrames

S_common = Dict(
    :size => (900,650),
    :background_color => :white,

    :legendfontsize=>13,

    :markerstrokewidth => 0,
    :grid => true,
    :gridstyle => :dash,
    :gridlinewidth =>2,
    :gridalpha=> 0.1,
)
S1_L = merge(S_common, Dict(
    :xlims => (0.4,2.9),
    :xticks=>(0.5:0.25:3),
    :ylims => (1, 6),
    :legend=>:topleft,
    :fname=>L"f_1(x) = ctg(x) + x^2"
))
S2_L = merge(S_common, Dict(
    :xlims => (-2.5,2.3),
    :xticks =>(-2.5:0.25:2),
    :ylims => (-10, 25),
    :legend=>:topright,
    :fname=>L"f_2(x) = x^5 - 3.2x^3 + 2.5x^2 - 7x + 1.5"
))
S1_H = merge(S_common, Dict(
    :xlims => (0.4,2.9),
    :xticks=>(0.5:0.25:3),
    :ylims => (-3.8, 7),
    :legend=>:topleft,
    :fname=>L"f_1(x) = ctg(x) + x^2"
))
S2_H = merge(S_common, Dict(
    :xlims => (-2.5, 2.2),
    :xticks =>(-2.5:0.25:2),
    :ylims => (-18,40),
    :legend=>:top,
    :fname=>L"f_2(x) = x^5 - 3.2x^3 + 2.5x^2 - 7x + 1.5"
))

nNodes = [4, 5, 6, 7, 10, 15, 30, 50]
fname1 = "CSVs/f1_$(nNodes[end]).csv"
df1 = CSV.read(fname1, DataFrame)
fname2 = "CSVs/f2_$(nNodes[end]).csv"
df2 = CSV.read(fname2, DataFrame)

function plotLagr(df, S)
    p = plot()

    # title!("Интерполяция  "*S[:fname][1:8]*"\$ полиномом Лагранжа",
    #     titlefontsize=S[:legendfontsize]+1
    # )

    plot!(
        df.x_exact, df.y_exact,
        label=S[:fname],
        color = :blue,
        line=:dash,
        linewidth = 1,
        markersize=3,
        ; S...
    )

    scatter!(
        df.nodesX, df.nodesY,
        label=latexstring("Узлы\\ функции:\\ $(count(!ismissing, df.nodesX))\\ шт."),
        color = :black,
        markeralpha=0.8,
        markersize=4,
        ; S...
    )

    plot!(
        df.x_exact, df.LagrangeY,
        label=L"Полином\ Лагранжа",
        color = :red,
        linewidth = 1.2,
        ; S...
    )

    plot!(twinx(),
        df.x_exact, df.errL,
        yaxis = (L"\varepsilon", (1e-16, 1e16), :log),
        xlims = S[:xlims],
        color=:black,
        alpha=0.35,
        label=false,
        yticks=((10. .^(-16:2:0))),
        ; S...
        
    )

    plot!([0],[0],
        color=:black,
        alpha=0.35,
        label=L"Профиль\ абсолютной\ ошибки\ \varepsilon"
    )
end

function plotHerm(df, S)
    plot()
    # title!("Интерполяция "*S[:fname][1:8]*"\$ сплайном Эрмита", titlefontsize=S[:legendfontsize]+1)
    plot!(
        df.x_exact, df.y_exact,
        label=S[:fname],
        color = :blue,
        line=:dash,
        linewidth = 1,
        markersize=3,
        ; S...
    )
    scatter!(
        df.nodesX, df.nodesY,
        label=latexstring("Узлы\\ функции:\\ $(count(!ismissing, df.nodesX))\\ шт."),
        color = :black,
        markeralpha=0.8,
        markersize=4,
        ; S...
    )
    plot!(
        df.x_exact, df.dfdx_exact,
        label=latexstring("f'_{$(S[:fname][4])}(x)"),
        color = :green,
        linewidth=8,
        alpha=0.25,
        markersize=3,
        ; S...
    )
    scatter!(
        df.nodesX, df.derNumY,
        label=latexstring("\\mathtt{f'}_{$(S[:fname][4])}(x)"),
        color = :green,
        markersize=6,
        alpha=.4,
        ; S...
    )
    plot!(
        df.nodesX, df.derNumY,
        label=false,
        line=:dash,
        color = :green,
        alpha=.4,
        ; S...
    )
    plot!(
        df.HermitX, df.HermitY,
        label=L"Сплайн\ Эрмита",
        color = :red,
        linewidth = 1.2,
        ; S...
    )
    
    plot!(twinx(),
        df.HermitX, df.errH,
        yaxis = (L"\varepsilon", (1e-16, 1e16), :log),
        color=:black,
        alpha=0.35,
        label=false,
        yticks=((10. .^(-16:2:0))),
        ; S...
    )
    plot!([0],[0],color=:black, alpha=0.35, label=L"Профиль\ абсолютной\ ошибки\ \varepsilon")
end

function plotsGen(nNodes, S1_L, S2_L, S1_H, S2_H)
    for n in nNodes
        filename1 = "CSVs/f1_$n.csv"
        filename2 = "CSVs/f2_$n.csv"
    
        df1 = CSV.read(filename1, DataFrame)
        df2 = CSV.read(filename2, DataFrame)
    
        l1 = plotLagr(df1, S1_L)
        l2 = plotLagr(df2, S2_L)
    
        h1 = plotHerm(df1, S1_H)
        h2 = plotHerm(df2, S2_H)
    
        savefig(l1, "plots/f1_Lagrange_$n.png")
        savefig(l2, "plots/f2_Lagrange_$n.png")
    
        savefig(h1, "plots/f1_Hermit_$n.png")
        savefig(h2, "plots/f2_Hermit_$n.png")
    end
end



# tex = "CSVs/f1_4.csv"

# dtex = CSV.read(tex, DataFrame)
# plotLagr(dtex, S1_L)
# plotHerm(dtex, S1_H)

# plotLagr(dtex, S2_L)
# plotHerm(dtex, S2_H)

# dtex[5:40:end, [5, 6, 9, 12, 13]][[1,3,4,5,6,7,8], :].errH  |> extrema|> print
# dtex.y_exact[5:40:end]

# dtex.derNumY[1:4] .- dtex.dfdx_exact[[1, 134, 267, 400]]

cheb = "CSVs/f2_Cheb_30.csv"
chebdf = CSV.read(cheb, DataFrame)

plotLagr(chebdf, S2_L)