# PIC8.PNG

using Plots, LaTeXStrings
using Plots: plot, plot!
using CSV
using DataFrames

S_common = Dict(
    :size => (1000, 680),
    :background_color => :white, :legendfontsize => 13, :markerstrokewidth => 0,
    :grid => true,
    :gridstyle => :dash,
    :gridlinewidth => 2,
    :gridalpha => 0.2,
)
S_err = merge(S_common, Dict(
    :yaxis => (L"\varepsilon", (1e-16, 1e10), :log),
    :xaxis => (L"N", (0, 82)),
    :legend => :topleft,
    :xticks => ((5:5:80)),
    :yticks => ((10.0 .^ (-16:2:10))),
    :linewidth => 2.5,
    # :alpha=>0.8,
))
apath = "~/Documents/Edu/4сем/Nummethods/Labs/3сем/0-2/tex/pics/"


function plotErr23(df1, df2, rand_dot_err_fname, S)
    p = plot()
    plot!(
        df1[!, "nNodes"], df1[!, "err_LagrUniform"],
        label=L"Трансцендентная\ -\ максимальная\ ошибка",
        color=:royalblue1,
        # line=:dash,
        ; S...
    )
    plot!(
        df2[!, "nNodes"], df2[!, "err_LagrUniform"],
        label=L"Полином\ -\ максимальная\ ошибка",
        color=:tomato,
        # line=:dash,
        ; S...
    )

    err_rand = CSV.read(rand_dot_err_fname, DataFrame)
    err_r_1 = err_rand.err_1
    err_r_2 = err_rand.err_2
    nodes = err_rand.nNodes
    plot!(
        nodes, err_r_1,
        label=L"Трансцендентная\ -\ ошибка\ в\ случайной\ точке",
        color=:royalblue1,
        line=:dash,
        ; S...
    )
    plot!(
        nodes, err_r_2,
        label=L"Полином\ -\ ошибка\ в\ случайной\ точке",
        color=:tomato,
        line=:dash,
        ; S...
    )

end


let
    maxerr_1 = Float64[]
    maxerr_2 = Float64[]
    for i=4:100
        e1 = CSV.read("CSVs/task2-3_f1_$i.csv", DataFrame).err_uniformLagr
        e2 = CSV.read("CSVs/task2-3_f2_$i.csv", DataFrame).err_uniformLagr
        push!(maxerr_1, e1 |> maximum)
        push!(maxerr_2, e2 |> maximum)
    end

    df1 = DataFrame("nNodes"=>4:100, "err_LagrUniform"=>maxerr_1)
    df2 = DataFrame("nNodes"=>4:100, "err_LagrUniform"=>maxerr_2)

    p8 = plotErr23(df1, df2, "CSVs/rdot_lagrange.csv", S_err)
    savefig(p8, apath * "pic8.png")
end