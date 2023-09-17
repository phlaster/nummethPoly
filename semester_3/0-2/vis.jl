using Plots, LaTeXStrings
using Plots: plot, plot!
using CSV
using DataFrames

function plotLagr(df, S)
    p = plot()
    plot!(
        df.x_uniformTab, df.y_uniformTab,
        label=S[:fname],
        color=:tomato,
        # line=:dash,
        linewidth=1.5,
        markersize=3,
        ; S...
    )

    scatter!(
        df.x_uniformNodes, df.y_uniformNodes,
        label=latexstring("Узлы\\ функции:\\ $(count(!ismissing, df.x_uniformNodes))\\ шт."),
        color=:black,
        markeralpha=0.8,
        markersize=4,
        ; S...
    )

    plot!(
        df.x_uniformTab, df.y_uniformLagr,
        label=L"Полином\ Лагранжа",
        color=:royalblue1,
        line=:dash,
        linewidth=1.2,
        ; S...
    )

    plot!(twinx(),
        df.x_uniformTab, df.err_uniformLagr,
        yaxis=(L"\varepsilon", (1e-16, 1e16), :log),
        xlims=S[:xlims],
        color=:black,
        alpha=0.35,
        label=false,
        yticks=((10.0 .^ (-16:2:0))),
        ; S...)

    plot!([0], [0],
        color=:black,
        alpha=0.35,
        label=L"Профиль\ абсолютной\ ошибки\ \varepsilon"
    )
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

function drawer23L(fname, N, S)
    fname = "CSVs/task2-3_$(fname)_$N.csv"
    df = CSV.read(fname, DataFrame)
    plotLagr(df, S)
end

function plotErr23(df1, df2, S)
    p = plot()
    plot!(
        df1[!, "nNodes"], df1[!, "err_LagrUniform"],
        label=L"Максимальная\ ошибка\ полинома\ Лагранжа\ для\ f_1(x)",
        color=:royalblue1,
        line=:dash,
        ; S...
    )
    plot!(
        df2[!, "nNodes"], df2[!, "err_LagrUniform"],
        label=L"Максимальная\ ошибка\ полинома\ Лагранжа\ для\ f_2(x)",
        color=:tomato,
        line=:dash,
        ; S...
    )
end

function plotHerm(df, S)
    plot()
    plot!(
        df.x_uniformTab, df.y_uniformTab,
        label=S[:fname],
        color=:tomato,
        line=:dash,
        linewidth=1.5,
        markersize=3,
        ; S...
    )
    scatter!(
        df.x_uniformNodes, df.y_uniformNodes,
        label=latexstring("Узлы\\ функции:\\ $(count(!ismissing, df.x_uniformNodes))\\ шт."),
        color=:black,
        markeralpha=0.8,
        markersize=4,
        ; S...
    )
    plot!(
        df.x_uniformTab, df.y_uniformTab_Der,
        label=latexstring("f'_{$(S[:fname][4])}(x) - точная\\ производная"),
        color=:mediumseagreen,
        linewidth=8,
        alpha=0.4,
        markersize=3,
        ; S...
    )
    scatter!(
        df.x_uniformNodes, df.y_uniformNodes_DerNum,
        label=latexstring("\\mathtt{f'}_{$(S[:fname][4])}(x) - метод\\ секущих"),
        color=:mediumseagreen,
        markersize=6,
        alpha=0.55,
        ; S...
    )
    plot!(
        df.x_uniformNodes, df.y_uniformNodes_DerNum,
        label=false,
        line=:dash,
        color=:mediumseagreen,
        alpha=0.55,
        ; S...
    )
    plot!(
        df.x_uniformTab, df.y_uniformHerm,
        label=L"Сплайн\ Эрмита",
        color=:royalblue1,
        linewidth=1.2,
        ; S...
    )

    plot!(twinx(),
        df.x_uniformTab, df.err_uniformHerm,
        yaxis=(L"\varepsilon", (1e-16, 1e16), :log),
        color=:black,
        alpha=0.35,
        label=false,
        yticks=((10.0 .^ (-16:2:0))),
        grid=S[:grid],
        gridstyle=S[:gridstyle],
        gridlinewidth=S[:gridlinewidth],
        gridalpha=S[:gridalpha],
        xlims=S[:xlims],)
    plot!([0], [0], color=:black, alpha=0.35, label=L"Профиль\ абсолютной\ ошибки\ \varepsilon")
end

function drawer23H(fname, N, S)
    fname = "CSVs/task2-3_$(fname)_$N.csv"
    df = CSV.read(fname, DataFrame)
    plotHerm(df, S)
end

function plotErr23andHerm(df1, df2, S)
    p = plot()
    plot!(
        df1[!, "nNodes"], df1[!, "err_LagrUniform"],
        label=L"L(f_1)",
        color=:royalblue1,
        line=:dash,
        ; S...
    )
    plot!(
        df2[!, "nNodes"], df2[!, "err_LagrUniform"],
        label=L"L(f_2)",
        color=:tomato,
        line=:dash,
        ; S...
    )
    plot!(
        df1[!, "nNodes"], df1[!, "err_Herm"],
        label=L"H(f_1)",
        color=:royalblue1,
        ; S...
    )

    plot!(
        df2[!, "nNodes"], df2[!, "err_Herm"],
        label=L"H(f_2)",
        color=:tomato,
        legendtitle=L"Средняя\ ошибка\ интерполяции:",
        ; S...
    )
end

function plot4L_cheb(df, S)
    plot(; S...)
    plot!(
        df.x_uniformTab, df.y_uniformTab,
        label=S[:fname],
        color=:tomato,
        line=:dash,
        linewidth=1.5,
        markersize=3,
    )

    scatter!(
        df.x_uniformNodes, df.y_uniformNodes,
        label=latexstring("Узлы\\ равномерной\\ решётки:\\ $(count(!ismissing, df.x_uniformNodes))\\ шт."),
        color=:black,
        # markeralpha=0.8,
        markersize=4,
    )


    scatter!(
        df.x_chebNodes, df.y_chebNodes,
        label=latexstring("Узлы\\ решётки\\ Чебышёва:\\ $(count(!ismissing, df.x_uniformNodes))\\ шт."),
        color=:saddlebrown,
        # markeralpha=0.8,
        marker=:xcross,
        markersize=8,
    )

    plot!(
        df.x_uniformTab, df.y_uniformLagr,
        label=latexstring("\$L("*S[:fname][2:4]*")\\ на\\ равномерной\\ решётке\$"),
        color=:royalblue,
        alpha=0.45,
        # line=:dash,
        linewidth=1.5,
        ; S...
    )

    plot!(
        df.x_uniformTab, df.y_chebLagr,
        label=latexstring("\$L("*S[:fname][2:4]*")\\ на\\ решётке\\ Чебышёва\$"),
        color=:royalblue,
        # line=:dash,
        linewidth=2.2,
        ; S...
    )

    plot!(twinx(),
        df.x_uniformTab, df.err_uniformLagr,
        yaxis=(L"\varepsilon", (1e-16, 1e16), :log),
        xlims=S[:xlims],
        color=:black,
        alpha=0.35,
        line=:dash,
        linewidth=1.2,
        label=false,
        yticks=((10.0 .^ (-16:2:0))),
        ; S...)
    plot!([0], [0],
        color=:black,
        alpha=0.35,
        line=:dash,
        label=L"Профиль\ \varepsilon_U\ равномерной\ решётки"
    )


    plot!(twinx(),
        df.x_uniformTab, df.err_chebLagr,
        yaxis=(L"\varepsilon", (1e-16, 1e16), :log),
        xlims=S[:xlims],
        color=:sandybrown,
        # alpha=0.65,
        # line=:dashdot,
        linewidth=1.2,
        label=false,
        yticks=((10.0 .^ (-16:2:0))),
        # ; S...
    )
    plot!([0], [0],
        color=:sandybrown,
        # alpha=0.65,
        # line=:dashdot,
        label=L"Профиль\ \varepsilon_T\ решётки\ Чебышёва"
    )
    
end

function drawer4L(fname, N, S)
    fname = "CSVs/task4_$(fname)_$N.csv"
    df = CSV.read(fname, DataFrame)
    plot4L_cheb(df, S)
end

function plotErr4(df1, df2, S)
    plot(; S...)
    plot!(
        df1[!, "nNodes"], df1[!, "err_LagrUniform"],
        label=L"L_U(f_1)\ -\ Трансцендентная\ ф., равномерная",
        color=:royalblue,
        line=:dash,
        ; S...
    )
    plot!(
        df2[!, "nNodes"], df2[!, "err_LagrUniform"],
        label=L"L_U(f_2)\ -\ Полиномиальная\ ф., равномерная",
        color=:tomato,
        line=:dash,
        ; S...
    )

    plot!(
        df1[!, "nNodes"], df1[!, "err_LagrCheb"],
        label=L"L_T(f_1)\ -\ Трансцендентная\ ф., Чебышёв",
        color=:royalblue,
        # line=:dash,
        ; S...
    )
    plot!(
        df2[!, "nNodes"], df2[!, "err_LagrCheb"],
        label=L"L_T(f_2)\ -\ Полиномиальная\ ф., Чебышёв",
        legendtitle=L"Интерполяция\ ф.\ полиномом\ Лагранжа\ на\ решётках",
        color=:tomato,
        # line=:dash,
        ; S...
    )
end

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

function plotNoizeProgression(df1, df2)
    plot(
        legend=:topleft,
        legendcolumns=2,
        xaxis=(:log, L"\mathbf{K}", (7e-7,0.4)),
        xticks=(1e1 .^ [-6:0;]),

        yaxis=((1e-6, 1e3), :log, L"\langle\varepsilon\rangle"),
        yticks=(1e1 .^ [-8:5;]),
        
        
        size = (1100, 800),
        background_color = :white, legendfontsize = 15, markerstrokewidth = 0,
        grid = true,
        gridlinewidth = 1.3,
        gridalpha = 0.25,
        minorgrid = true,
        minorgridalpha = 0.3,
    )

    plot!(
        df1.dev,
        [df1.err_LagrUniform_Noize df2.err_LagrUniform_Noize df1.err_Herm_Noize df2.err_Herm_Noize df1.err_LagrCheb_Noize df2.err_LagrCheb_Noize],
        label=[L"L_U(f_1)" L"L_U(f_2)" L"H_U(f_1)" L"H_U(f_2)" L"L_T(f_1)" L"L_T(f_2)"],
        ls=[:solid :solid :solid :solid :dot :dot],
        color=[:royalblue3 :red3 :steelblue3 :tomato2 :skyblue :lightsalmon],
        linewidth=4,
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
S1_H = merge(S_common, Dict(
    :xlims => (0.4, 2.9),
    :xticks => (0.5:0.25:3),
    :ylims => (-2, 6),
    :yticks => (-2:0.5:6),
    :legend => :topleft,
    :fname => L"f_1(x) = ctg(x) + x^2"
))
S2_H = merge(S_common, Dict(
    :xlims => (-2.5, 2.2),
    :xticks => (-2.5:0.25:2),
    :ylims => (-18, 40),
    :legend => :top,
    :fname => L"f_2(x) = x^5 - 3.2x^3 + 2.5x^2 - 7x + 1.5"
))
S_err = merge(S_common, Dict(
    :yaxis => (L"\varepsilon", (1e-16, 1e10), :log),
    :xaxis => (L"N", (0, 82)),
    :legend => :topleft,
    :xticks => ((5:5:80)),
    :yticks => ((10.0 .^ (-16:2:10))),
    :linewidth => 2.5,
    # :alpha=>0.8,
))

if abspath(PROGRAM_FILE) == @__FILE__
    Lagr1 = CSV.read("CSVs/task2-3_f1_4.csv", DataFrame)
    Lagr1[[1, 2, 3, 399, 400], [4, 5, 7, 9]]
    apath = "/home/alex/Documents/Edu/4сем/Nummethods/Labs/Отчёты/0-2/pics/"
    Lagr1.err_uniformLagr |> maximum |> println
    
    p1 = drawer23L("f1", 4, S1_L)
    savefig(p1, apath * "pic1.png")
    p2 = drawer23L("f1", 9, S1_L)
    savefig(p2, apath * "pic2.png")
    p3 = drawer23L("f1", 30, S1_L)
    savefig(p3, apath * "pic3.png")
    p4 = drawer23L("f1", 50, S1_L)
    savefig(p4, apath * "pic4.png")
    errsdf1 = CSV.read("CSVs/task2-3_f1_summ_4-100_.csv", DataFrame)
    errsdf1[[1, 7, 27, 47], :]
    p5 = drawer23L("f2", 4, S2_L)
    savefig(p5, apath * "pic5.png")
    p6 = drawer23L("f2", 6, S2_L)
    savefig(p6, apath * "pic6.png")
    p7 = drawer23L("f2", 30, S2_L)
    savefig(p7, apath * "pic7.png")
    errsdf2 = CSV.read("CSVs/task2-3_f2_summ_4-100_.csv", DataFrame)
    errsdf2[[1, 3, 27], :]

    df1 = CSV.read("CSVs/task2-3_f1_summ_4-100_.csv", DataFrame)
    df2 = CSV.read("CSVs/task2-3_f2_summ_4-100_.csv", DataFrame)
    
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

        p8 = plotErr23(df1, df2, S_err)
        savefig(p8, apath * "pic8.png")
    end

    p9 = drawer23H("f1", 4, S1_H)
    savefig(p9, apath * "pic9.png")
    p10 = drawer23H("f1", 9, S1_H)
    savefig(p10, apath * "pic10.png")
    p11 = drawer23H("f1", 30, S1_H)
    savefig(p11, apath * "pic11.png")
    p12 = drawer23H("f1", 50, S1_H)
    savefig(p12, apath * "pic12.png")
    errsdf1 = CSV.read("CSVs/task2-3_f1_summ_4-100_.csv", DataFrame)
    errsdf1[[1, 7, 27, 47], :]
    Herm1 = CSV.read("CSVs/task2-3_f1_4.csv", DataFrame)
    Herm1[[1, 2, 3, 399, 400], [4, 5, 8, 10]]
    errsdf1[[1, 7, 27, 47], :]

    p13 = drawer23H("f2", 4, S2_H)
    savefig(p13, apath * "pic13.png")
    p14 = drawer23H("f2", 7, S2_H)
    savefig(p14, apath * "pic14.png")
    p15 = drawer23H("f2", 15, S2_H)
    savefig(p15, apath * "pic15.png")

    df1 = CSV.read("CSVs/task2-3_f1_summ_4-100_.csv", DataFrame)
    df2 = CSV.read("CSVs/task2-3_f2_summ_4-100_.csv", DataFrame)
    p16 = plotErr23andHerm(df1, df2, S_err)
    savefig(p16, apath * "pic16.png")

    p18 = drawer4L("f1", 4, S1_L)
    savefig(p18, apath * "pic18.png")
    p19 = drawer4L("f1", 5, S1_L)
    savefig(p19, apath * "pic19.png")
    p20 = drawer4L("f1", 13, S1_L)
    savefig(p20, apath * "pic20.png")
    p21 = drawer4L("f1", 35, S1_L)
    savefig(p21, apath * "pic21.png")

    p22 = drawer4L("f2", 4, S2_L)
    savefig(p22, apath * "pic22.png")
    p23 = drawer4L("f2", 5, S2_L)
    savefig(p23, apath * "pic23.png")
    p24 = drawer4L("f2", 6, S2_L)
    savefig(p24, apath * "pic24.png")
    p25 = drawer4L("f2", 35, S2_L)
    savefig(p25, apath * "pic25.png")

    df1 = CSV.read("CSVs/task4_f1_summ_4-100_.csv", DataFrame)
    df2 = CSV.read("CSVs/task4_f2_summ_4-100_.csv", DataFrame)
    p26 = plotErr4(df1, df2, S_err)
    savefig(p26, apath * "pic26.png")

    p27 = plotNoize("f1_10_0.100000", S1_L)
    savefig(p27, apath * "pic27.png")

    p28 = plotNoize("f2_25_0.200000", S2_L)
    savefig(p28, apath * "pic28.png")

    df1 = CSV.read("CSVs/task5_f1_prog_20_0.200000_.csv", DataFrame)
    df2 = CSV.read("CSVs/task5_f2_prog_20_0.200000_.csv", DataFrame)

    p29 = plotNoizeProgression(df1, df2)
    savefig(p29, apath * "pic29.png")
end