using CSV, DataFrames


maxerr_1 = Float64[]
maxerr_2 = Float64[]

for i=4:100
    e1 = CSV.read("CSVs/task2-3_f1_$i.csv", DataFrame).err_uniformHerm
    e2 = CSV.read("CSVs/task2-3_f2_$i.csv", DataFrame).err_uniformHerm
    push!(maxerr_1, e1 |> maximum)
    push!(maxerr_2, e2 |> maximum)
end

maxerr_1[50-3]
maxerr_2[15-3]

function get_max_errs(taskf="task2-3_f1", col=:err_uniformLagr; lims=4:81)
    maxerrs = Float64[]

    for i=lims
        e = CSV.read("CSVs/$(taskf)_$i.csv", DataFrame)[!, col]
        push!(maxerrs, e |> maximum)
    end
    
    maxerrs
end

get_max_errs()

get_max_errs.("task4_f2", [:err_uniformLagr, :err_chebLagr]) .|> (x->x[6-3])