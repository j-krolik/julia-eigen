using BenchmarkTools
using Plots

mutable struct TestedFunction
    fun::Function
    name::String
    result
    TestedFunction(fun, name) = new(fun, name, nothing)
end


# create tested_functions list
tested_functions = []
push!(tested_functions, TestedFunction(eigen,"default"))
push!(tested_functions, TestedFunction(eigen_power,"power"))
push!(tested_functions, TestedFunction(eigen_qr,"qr"))
push!(tested_functions, TestedFunction(eigen_householder,"house"))

# clear tested_functions
for i in 1:size(tested_functions)[1]
    pop!(tested_functions)
end

foo = [12 -51 4; 6 167 -68; -4 24 -41] * 1.0

function bench(a=foo)
    global fun
    for tested_fun in tested_functions
        fun = tested_fun.fun
        b = @benchmarkable (fun($a)) seconds=1.0
        tune!(b)
        tested_fun.result = (run(b))
    end
    print_bench()
end

function print_bench()
    for tested_fun in tested_functions
        println(tested_fun.result)
    end
end

##################
##################


t_bench_group = BenchmarkGroup()
t_functions = (eigen_power, eigen_qr, eigen_householder)
#t_functions = (eigen_power, eigen_householder)
t_matrices = (3,3), (5,5) #, (20,20), (50,50), (100,100))
t_iteration_no = 3

t_functions_string = Vector{String}(undef, size(t_functions,1))
for f_no = 1:size(t_functions,1)
    t_functions_string[f_no] = string(t_functions[f_no])
end
t_matrices_string = Vector(undef, size(t_matrices,1))
for matix_no = 1:size(t_matrices,1)
    t_matrices_string[matix_no] = string(size(t_matrices[matix_no],1))
end


t_A = Array{Any,2}(undef, t_iteration_no, size(t_matrices,1))
for i = 1:t_iteration_no, A_no = 1:size(t_matrices,1)
    t_A[i,A_no] = rand(t_matrices[A_no][1],t_matrices[A_no][2])
end


for f in t_functions
    t_bench_group[string(f)] = BenchmarkGroup([f])
    for A_no = 1:size(t_matrices,1)
        t_bench_group[string(f), string(size(A))] = BenchmarkGroup([f,A])
        for i = 1:t_iteration_no
            A = t_A[i,A_no]
            t_bench_group[string(f), string(size(A)), i] = @benchmarkable $(f)($A)
        end
    end
end

#t_bench_group[["eigen_power","(3, 3)"]]
#t_bench_group[[string(t_functions[1]),string(size(t_matrices[1]))]]

# tune!(t_bench_group)
t_results_raw = run(t_bench_group, verbose = true, seconds = 1)

# t_results[string(t_functions[1])][string(size(t_matrices[1]))]
# t_results[[string(t_functions[1]), string(size(t_matrices[1]))]]
# median(t_results_raw[string(t_functions[1])][string(size(t_matrices[1]))])

t_result_time = Matrix{Float64}(undef,size(t_matrices,1),size(t_functions,1))
t_result_memory = Matrix{Int64}(undef,size(t_matrices,1),size(t_functions,1))
for f_no = 1:size(t_functions,1)
    for A_no = 1:size(t_matrices,1)
        time_sum   = 0
        memory_sum = 0
        for i = 1:t_iteration_no
            median_result = median(t_results_raw[string(t_functions[f_no])][string(size(t_matrices[A_no]))])
            time_sum   += median_result.time
            memory_sum += median_result.memory
        end
        t_result_time[A_no,f_no]   = time_sum/t_iteration_no
        t_result_memory[A_no,f_no] = memory_sum/t_iteration_no
    end
end

plot(t_matrices_string, t_result_time/1000,
#plot(t_matrices_string, t_result_time[:,[1,3]],
    yaxis = (:log10),
    title = "Mediana czasu w zależności wyznacznia wektrów od rozmiaru macierzy",
    titlefontsize = 10,
    ylabel = "Czas, ms",
    xlabel = "Rozmar macierzy",
    legend = :topleft,
    label = ["Power" "QR" "Householder"])
    #label = ["Power" "Householder"])

plot(t_matrices_string, t_result_memory/1024,
#plot(t_matrices_string, t_result_time[:,[1,3]],
    yaxis = (:log10),
    title = "Użycie pamięci w zależności od rozmiaru macierzy",
    titlefontsize = 10,
    ylabel = "Użycie pamięci, KiB",
    xlabel = "Rozmar macierzy",
    legend = :topleft,
    label = ["Power" "QR" "Householder"])
    #label = ["Power" "Householder"])


plot(["a","c","b"], 1:3)
