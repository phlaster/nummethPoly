using CairoMakie, CSV, ColorSchemes

f1 = CSV.File("rect_trapez_compare.csv");
f_ru = CSV.File("trapez_runge.csv");
f_ad = CSV.File("trapez_adaptive.csv");

let
    f = Figure(resolution = (1000, 800))
    ax = CairoMakie.Axis(f[1, 1],
        backgroundcolor = :gray90,
        xlabel = "Шаги", ylabel = "абс. ошибка",
        xminorgridvisible = true,
        yminorgridvisible = true,
        xscale=log10,
        yscale=log10,
	)
    lines!(
        f1.n_steps, f1.err_rect,
        color = :red,
        label = "прямоугольный (центральная точка)",
        linestyle = :dot
    )
    scatter!(
        f1.n_steps, f1.err_trapez,
        color = :blue,
        label = "трапеции",
        
    )
    scatter!(
		f1.n_steps, replace(f1.delta_err, 0=>missing),
        color = :gray60,
        label = "различие",
	)
    
    axislegend()
    current_figure()
end

let
    f = Figure(resolution = (1000, 800))
    ax = CairoMakie.Axis(f[1, 1],
        backgroundcolor = :gray90,
        xlabel = "Шаги", ylabel = "абс. ошибка",
        xminorgridvisible = true,
        yminorgridvisible = true,
        xscale=log10,
        yscale=log10,
	)
    # lines!(
    #     f1.n_steps, f1.err_rect,
    #     color = :red,
    #     label = "прямоугольный (центральная точка)"
    # )
    scatter!(
        f2.n_steps, f2.err,
        color = :royalblue,
        label = "трапеции (Рунге)",
        markersize = 25,
        alpha=.7
    )
    scatter!(
        f1.n_steps, f1.err_trapez,
        color = :blue,
        label = "трапеции",
        markersize = 10,
    )
    scatter!(
		f3.n_steps, f3.err,
        color = :cyan3,
        label = "трапеции адаптивный",
	)
    
    axislegend()
    current_figure()
end

let
    f = Figure(resolution = (1000, 800))
    ax = CairoMakie.Axis(f[1, 1],
        backgroundcolor = :gray90,
        ylabel = "Шаги", xlabel = "точность",
        # xminorgridvisible = true,
        # yminorgridvisible = true,
        xscale=log10,
        yscale=log10,
	)
    scatter!(
        f2.eps, f2.n_steps,
        color = :royalblue,
        label = "трапеции (Рунге)",
        markersize = 25,
        alpha=.7
    )
    scatter!(
        f3.eps, f3.n_steps,
        color = :cyan3,
        label = "трапеции адаптивный"
    )
    # ylims!(1e-15, 1)    
    axislegend()
    current_figure()
end

let
    f = Figure(resolution = (1000, 800))
    ax = CairoMakie.Axis(f[1, 1],
        backgroundcolor = :gray90,
        xlabel = "точность", ylabel = "абс. ошибка",
        xminorgridvisible = true,
        yminorgridvisible = true,
        xscale=log10,
        yscale=log10,
	)
    scatter!(
        f2.eps, f2.err,
        color = :royalblue,
        label = "трапеции (Рунге)",
        markersize = 25,
        alpha=.7
    )
    scatter!(
        f3.eps, f3.err,
        color = :cyan3,
        label = "трапеции адаптивный"
    )
    xlims!(1e-13, 10)
    ylims!(1e-13, 10) 
    axislegend()
    current_figure()
end

let 
    f = Figure(resolution = (1200, 800))
    ax1 = CairoMakie.Axis(f[1, 1],
        backgroundcolor = :gray90,
        xlabel = "ε", ylabel = "N",
        # xminorgridvisible = true,
        # yminorgridvisible = true,
        xscale=log10,
        yscale=log10,
	)
    scatterlines!(
        ax1, f_ru.eps, f_ru.n_steps,
        label = "Трапеции, критерий Рунге"
    )
    scatterlines!(
        ax1, f_ad.eps, f_ad.n_steps,
        label = "Адаптивные трапеции",
        color =:crimson
    )
    axislegend()


    ax2 = CairoMakie.Axis(f[2, 1],
        backgroundcolor = :gray90,
        xlabel = "ε", ylabel = "err",
        # xminorgridvisible = true,
        # yminorgridvisible = true,
        xscale=log10,
        yscale=log10,
	)
    scatterlines!(
        ax2, f_ru.eps, f_ru.err,
        label = "Трапеции, критерий Рунге"
    )
    scatterlines!(
        ax2, f_ad.eps, f_ad.err,
        label = "Адаптивные трапеции",
        color =:crimson
    )
    # axislegend(position=:lt)
    current_figure()
end