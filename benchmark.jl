using BenchmarkTools
using Plots

bench_group = BenchmarkGroup()
iteration_no = 21
functions = [eigen_power, eigen_qr, eigen_householder]
matrices = (3,3),(5,5),(10,10),(20,20),(40,40),(60,60),(80,80),(100,100)
matrices_string = Vector(undef, size(matrices,1))
for matix_no = 1:size(matrices,1)
    matrices_string[matix_no] = string(matrices[matix_no])
end

# generowanie macierzy, które będą podawane benchmarkowi
A = Array{Any,2}(undef, iteration_no, size(matrices,1))
for i = 1:iteration_no, A_no = 1:size(matrices,1)
    A[i,A_no] = rand(matrices[A_no][1],matrices[A_no][2])
end

# definiowanie procesu benchmarku
for f in functions
    bench_group[string(f)] = BenchmarkGroup([f])
    for A_no = 1:size(matrices,1)
        bench_group[string(f), string(size(A))] = BenchmarkGroup([f,A])
        for i = 1:iteration_no
            A = A[i,A_no]
            bench_group[string(f), string(size(A)), i] = @benchmarkable $(f)($A)
        end
    end
end

# uruchomienie benchmarka
results_raw = run(bench_group, verbose = true, seconds = 1)

# ekstracja danych z postaci surowej
resultime = Matrix{Float64}(undef,size(matrices,1),size(functions,1)+1)
resulmemory = Matrix{Float64}(undef,size(matrices,1),size(functions,1)+1)
for f_no = 1:size(functions,1)
    for A_no = 1:size(matrices,1)
        time_sum   = Float64[]
        memory_sum = Float64[]
        for i = 1:iteration_no
            if f_no < 4
                median_result = median(results_raw[string(functions[f_no]), string(matrices[A_no]), i])
            else
                median_result = median(results_raw_small[string(functions[f_no-3]), string(matrices[A_no]), i])
            end
            push!(time_sum, median_result.time)
            push!(memory_sum, median_result.memory)
        end
        resultime[A_no,f_no]   = sum(time_sum)/iteration_no
        resulmemory[A_no,f_no] = sum(memory_sum)/iteration_no
    end
end

plot(matrices_string, resultime/1000,
    yaxis = (:log10),
    title = "Mediana czasu w zależności wyznacznia wektrów od rozmiaru macierzy",
    titlefontsize = 10,
    ylabel = "Czas, ms",
    xlabel = "Rozmar macierzy",
    legend = :topleft,
    label = ["Power Δ=10^-7" "QR Δ=10^-7" "Householder Δ=10^-7" "Power Δ=10^-6"])

plot(matrices_string, resulmemory/1024,
    yaxis = (:log10),
    title = "Użycie pamięci w zależności od rozmiaru macierzy",
    titlefontsize = 10,
    ylabel = "Użycie pamięci, KiB",
    xlabel = "Rozmar macierzy",
    legend = :topleft,
    label = ["Power Δ=10^-7" "QR Δ=10^-7" "Householder Δ=10^-7" "Power Δ=10^-6"])
