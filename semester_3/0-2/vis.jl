using Plots, ColorSchemes, LaTeXStrings
using CSV
using DataFrames

fname1 = "f1.csv"
df1 = CSV.read(fname1, DataFrame)
fname2 = "f2.csv"
df2 = CSV.read(fname2, DataFrame)

S1 = Dict(
    :size => (800,500),
    :xlims => (0.4,2.9),
    :ylims => (-5.3, 7),
    :background_color => :white,
    :legend=>:bottom,
    # :legendfontsize=>13,
    
    :grid => true,
    :gridalpha=> 0.7,

    :minorgrid => true,
    :minorgridalpha=> 0.2,

    
    # :markeralpha=>0.5,
    :markerstrokewidth => 0,

    # :linewidth => 3,
    # :linealpha => 0.7,
    
    # :palette => ColorSchemes.berlin10[1:3:end],
    :fname=>L"f_1(x) = ctg(x) + x^2"
)

S2 = Dict(
    :size => (800,500),
    :xlims => (-2.5, 2.2),
    :ylims => (-18,55),
    :background_color => :white,
    :legend=>:top,
    # :legendfontsize=>13,
    
    :grid => true,
    :gridalpha=> 0.7,

    :minorgrid => true,
    :minorgridalpha=> 0.2,

    :markerstrokewidth => 0,

    :fname=>L"f_2(x) = x^5 - 3.2x^3 + 2.5x^2 - 7x + 1.5"
)

function plotall(df, S)
    p = scatter(
        df.nodesX, df.nodesY,
        label="Узлы ($(count(!ismissing, df.nodesX))шт.)",
        color = :black,
        markeralpha=0.4,
        markersize=7.7,

        ; S...
    )
    title!(S[:fname])
    plot!(
        df.x_exact, df.y_exact,
        label=L"f(x)"*" точно",
        color = :black,
        alpha = 0.45,
        linewidth = 8,
        markersize=3,
        ; S...
    )
    plot!(
        df.x_exact, df.dfdx_exact,
        label=L"f'(x)"*" точно",
        color = :green,
        linewidth=8,
        alpha=0.25,
        markersize=3,
        ; S...
    )
    scatter!(
        df.nodesX, df.derNumY,
        label=L"f'(x)"*" численно по узлам",
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
        df.LagrangeX, df.LagrangeY,
        label="Полином Лагранжа",
        color = :blue,
        line=:dash,
        linewidth = 4,
        ; S...
    )
    plot!(
        df.HermitX, df.HermitY,
        label="Сплайн",
        color=:red,
        line=:dashdot,
        linewidth=2,
    )
    
    plot!(twinx(),
        df.LagrangeX, df.errL,
        yaxis = ((1e-16, maximum(df.errH)*1e12), :log),
        minorgrid=true,
        yticks=((10. .^(-16:2:0))),
    )

    return p
end



function plotHerm(df, S)

end


p1 = plotall(df1, S1)

p2 = plotall(df2, S2)



function plotLagr(df, S)
    plot()
    title!("Интерполяция  "*S[:fname][1:8]*"\$ полиномом Лагранжа", titlefontsize=S[:legendfontsize]+1)
    plot!(
        df.x_exact, df.y_exact,
        label=S[:fname],
        color = :blue,
        # alpha = 0.45,
        line=:dash,
        linewidth = 1,
        markersize=3,
        ; S...
    )
    scatter!(
        df.nodesX, df.nodesY,
        label="Узлы функции: $(count(!ismissing, df.nodesX)) шт.",
        color = :black,
        markeralpha=0.8,
        markersize=4,
        ; S...
    )
    plot!(
        df.LagrangeX, df.LagrangeY,
        label="Полином Лагранжа",
        color = :red,
        line=:dashdot,
        linewidth = 1,
        ; S...
    )
    
    plot!(twinx(),
        df.LagrangeX, df.errL,
        yaxis = (L"ε", (1e-16, maximum(df.errL)*1e12), :log),
        xlims = S[:xaxis][2],
        color=:grey60,
        minorgrid=true,
        legend=false,
        yticks=((10. .^(-16:2:0))),
    )
    plot!([0],[0],color=:grey60, label=L"ε")
end
S1_L = Dict(
    :size => (700,500),
    :xaxis=> (L"x", (0.4,2.9)),
    :yaxis=> (L"f_1(x)", (1, 6)),
    # :xlims => (0.4,2.9),
    # :ylims => (0, 6),
    :background_color => :white,
    :legend=>:topleft,
    :legendfontsize=>11,
    
    :grid => true,
    :gridalpha=> 0.7,

    :minorgrid => true,
    :minorgridalpha=> 0.2,

    
    # :markeralpha=>0.5,
    :markerstrokewidth => 0,

    # :linewidth => 3,
    # :linealpha => 0.7,
    
    # :palette => ColorSchemes.berlin10[1:3:end],
    :fname=>L"f_1(x) = ctg(x) + x^2"
)
plotLagr(df1, S1_L)

S2_L = Dict(
    :size => (700,500),
    :xaxis=> (L"x", (-2.5,2.3)),
    :yaxis=> (L"f_2(x)", (-10, 25)),
    # :xlims => (0.4,2.9),
    # :ylims => (0, 6),
    :background_color => :white,
    :legend=>:topright,
    :legendfontsize=>11,
    
    :grid => true,
    :gridalpha=> 0.7,

    :minorgrid => true,
    :minorgridalpha=> 0.2,

    
    # :markeralpha=>0.5,
    :markerstrokewidth => 0,

    # :linewidth => 3,
    # :linealpha => 0.7,
    
    # :palette => ColorSchemes.berlin10[1:3:end],
    :fname=>L"f_2(x) = x^5 - 3.2x^3 + 2.5x^2 - 7x + 1.5"
)
plotLagr(df2, S2_L)


function plotHerm(df, S)
    plot()
    title!("Интерполяция "*S[:fname][1:8]*"\$ сплайном Эрмита", titlefontsize=S[:legendfontsize]+1)
    plot!(
        df.x_exact, df.y_exact,
        label=S[:fname],
        color = :blue,
        # alpha = 0.45,
        line=:dash,
        linewidth = 1,
        markersize=3,
        ; S...
    )
    scatter!(
        df.nodesX, df.nodesY,
        label="Узлы функции: $(count(!ismissing, df.nodesX)) шт.",
        color = :black,
        markeralpha=0.8,
        markersize=4,
        ; S...
    )
    plot!(
        df.x_exact, df.dfdx_exact,
        label=L"f'(x)"*" точно",
        color = :green,
        linewidth=8,
        alpha=0.25,
        markersize=3,
        ; S...
    )
    scatter!(
        df.nodesX, df.derNumY,
        label=L"f'(x)"*" численно по узлам",
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
        label="Сплайн Эрмита",
        color = :red,
        line=:dash,
        linewidth = 1,
        ; S...
    )
    
    plot!(twinx(),
        df.HermitX, df.errH,
        yaxis = (L"ε", (1e-16, maximum(df.errH)*1e12), :log),
        xlims = S[:xaxis][2],
        color=:grey60,
        minorgrid=true,
        legend=false,
        yticks=((10. .^(-16:2:0))),
    )
    plot!([0],[0],color=:grey60, label=L"ε")
end

S1_H = Dict(
    :size => (700,500),
    :xaxis=> (L"x", (0.4,2.9)),
    :yaxis=> (L"f_1(x)", (-4, 8)),
    :xlims => (0.4,2.9),
    :ylims => (-5.3, 7),
    :background_color => :white,
    :legend=>:topleft,
    :legendfontsize=>11,
    
    :grid => true,
    :gridalpha=> 0.7,

    :minorgrid => true,
    :minorgridalpha=> 0.2,

    
    # :markeralpha=>0.5,
    :markerstrokewidth => 0,

    # :linewidth => 3,
    # :linealpha => 0.7,
    
    # :palette => ColorSchemes.berlin10[1:3:end],
    :fname=>L"f_1(x) = ctg(x) + x^2"
)
plotHerm(df1, S1_H)

S2_H = Dict(
    :size => (700,500),
    :xaxis=> (L"x", (-2.5,2.3)),
    :yaxis=> (L"f_2(x)", (-20, 50)),
    :xlims => (-2.5, 2.2),
    :ylims => (-18,55),
    :background_color => :white,
    :legend=>:top,
    :legendfontsize=>8,
    
    :grid => true,
    :gridalpha=> 0.7,

    :minorgrid => true,
    :minorgridalpha=> 0.2,

    
    # :markeralpha=>0.5,
    :markerstrokewidth => 0,

    # :linewidth => 3,
    # :linealpha => 0.7,
    
    # :palette => ColorSchemes.berlin10[1:3:end],
    :fname=>L"f_2(x) = x^5 - 3.2x^3 + 2.5x^2 - 7x + 1.5"
)
plotHerm(df2, S2_H)
