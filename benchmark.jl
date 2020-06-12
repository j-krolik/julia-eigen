using BenchmarkTools

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


t_bench_group = BenchmarkGroup()
t_functions = (eigen_power, eigen_qr, eigen_householder)
t_matrices = (rand(3,3), rand(5,5), rand(20,20) )#, rand(50,50), rand(100,100))

for f in t_functions
    t_bench_group[string(f)] = BenchmarkGroup([f])
    for A in t_matrices
        t_bench_group[string(f)][string(size(A))] = @benchmarkable $(f)($A)
    end
end
#t_bench_group[["eigen_power","(3, 3)"]]
#t_bench_group[[string(t_functions[1]),string(size(t_matrices[1]))]]

# tune!(t_bench_group)

t_results = run(t_bench_group, verbose = true, seconds = 1)


# t_results[string(t_functions[1])][string(size(t_matrices[1]))]
# t_results[[string(t_functions[1]), string(size(t_matrices[1]))]]
median(t_results[string(t_functions[1])][string(size(t_matrices[1]))])
