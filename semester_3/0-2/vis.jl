using Plots, ColorSchemes
using CSV
using DataFrames

fname1 = "f1.csv"
df1 = CSV.read(fname1, DataFrame)

S1 = Dict(
    :size => (800,400),
    :xlims => (0.2,3),
    :ylims => (-5, 8),
    :background_color => :grey90,
    
    :grid => true,
    :gridalpha=> 0.7,

    :minorgrid => true,
    :minorgridalpha=> 0.2,

    
    :markeralpha=>0.5,
    :markerstrokewidth => 0,

    :linewidth => 3,
    :linealpha => 0.7,
    
    # :palette => ColorSchemes.berlin10[1:3:end],
)

nodes = scatter(
    df1.nodesX, df1.nodesY,
    label="f1 (узлы)",
    color = :black,
    markersize=7,
    ; S1...
)

tabulate = scatter!(
    df1.x_exact, df1.y_exact,
    label="f1",
    color = :red,
    markersize=3,
    ; S1...
)
deriv = plot!(
    df1.x_exact, df1.dfdx_exact,
    label="Производная",
    color = :blue,
    ; S1...
)

lagrange = plot!(
    df1.LagrangeX, df1.LagrangeY,
    label="Полином Лагранжа",
    color = :green
    ; S1...
)


plot!(df1.HermitX, df1.HermitY)

fname2 = "f2.csv"
df2 = CSV.read(fname2, DataFrame)
scatter(df2.x_exact, df2.y_exact, label="f2")
scatter(df2.x_exact, df2.y_exact)
plot(df2.x_exact, df2.y_exact)
plot!(df2.x_exact, df2.dfdx_exact)
plot!(df2.x_exact, df2.Lagrange)
plot!(df2.HermitX, df2.HermitY)