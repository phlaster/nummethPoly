using CairoMakie, CSV, ColorSchemes

f1 = CSV.File("baza_grid=0.1.csv")
f2 = CSV.File("baza_grid=0.05.csv")
f3 = CSV.File("baza_h_i=1e-1_1e-8.csv")

baza_picture(file) = begin
    f = Figure(resolution = (1000, 800))
    ax1 = CairoMakie.Axis(f[1, 1],
        backgroundcolor = :gray90,
        xlabel = "x", ylabel = "y",
        xminorgridvisible = true,
        yminorgridvisible = true,
	)
    lines!(
        file.x_i, file.y_i,
        color = :red,
        label = "точное значение",
        )
    scatter!(
        file.x_i, file.y_appr_i,
        color = :royalblue,
        label = "метод Рунге-Кутты",
        markersize=8
    )
    
    ax2 = CairoMakie.Axis(f[2, 1],
    backgroundcolor = :gray90,
    xlabel = "x", ylabel = "eps",
    xminorgridvisible = true,
    yminorgridvisible = true,
    )
    scatter!(
        file.x_i, file.err,
        color = :gray50,
        label = "ошибка",
        markersize=8
    )
    ylims!(0, 0.3)  
    Legend(f[1, 2], ax1)
    Legend(f[2, 2], ax2)
    f
end

baza_picture(f1)
baza_picture(f2)
baza(CSV.File("test.csv"))

let
    f = Figure(resolution = (1000, 700))
    ax = CairoMakie.Axis(f[1, 1],
        backgroundcolor = :gray90,
        xlabel = "h_i", ylabel = "err",
        xscale=log10,
        yscale=log10,
	)
    scatter!(
        f3.h_i, f3.err,
        color = :royalblue,
        label = "метод Рунге-Кутты",
        markersize=8
    )
    xlims!(1e-9, 1e0)
    lines!(
        f3.h_i, f3.h_i .^2,
        color = :red,
        label = "h^2",
    )
    Legend(f[1, 2], ax)
    # Legend(f[2, 2], ax2)
    f
end