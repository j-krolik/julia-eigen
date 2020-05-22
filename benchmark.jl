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

function bench!(tested_function_list)
    global fun
    for tested_fun in tested_function_list
        fun = tested_fun.fun
        b = @benchmarkable (fun($foo))
        tune!(b)
        #tested_fun.result = 1
        tested_fun.result = median(run(b))
    end
end

function print_bench()
    for tested_fun in tested_functions
        println(tested_fun.result)
    end
end
