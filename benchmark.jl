using BenchmarkTools

mutable struct TestedFunction
    fun::Function
    name::String
    result
    TestedFunction(fun, name) = new(fun, name, nothing)
end

tested_functions = []
push!(tested_functions, TestedFunction(eigen,"default"))
push!(tested_functions, TestedFunction(start_power,"power"))
push!(tested_functions, TestedFunction(start_qr,"qr"))
push!(tested_functions, TestedFunction(start_house,"house"))

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
