using CairoMakie, CSV, ColorSchemes

f1 = CSV.File("rect_trapez_compare.csv");
f2 = CSV.File("trapez_runge.csv");
f3 = CSV.File("trapez_adaptive.csv");

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
        xlabel = "Шаги", ylabel = "абс. ошибка",
        xminorgridvisible = true,
        yminorgridvisible = true,
        xscale=log10,
        yscale=log10,
	)
    scatter!(
        f2.n_steps, f2.err,
        color = :royalblue,
        label = "трапеции (Рунге)",
        markersize = 25,
        alpha=.7
    )
    scatter!(
        f3.n_steps, f3.err,
        color = :cyan3,
        label = "трапеции адаптивный"
    )
    ylims!(1e-15, 1)    
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

