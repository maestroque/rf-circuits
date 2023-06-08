using Revise
using Plots
using Statistics
using Evolutionary
using Distributions
using LaTeXStrings
includet("transmissionLines.jl")

function stubCircuitReflection(lengthsVector, f)
    d1 = lengthsVector[1]
    l1 = lengthsVector[2]
    d2 = lengthsVector[3]
    l2 = lengthsVector[4]
    d3 = lengthsVector[5]
    l3 = lengthsVector[6]

    zLoad = 100
    z0 = 50
    f0 = 1e9
    lengthToLambda = l -> l * f * f0 / 300000000  

    zIn1 = zIn(z0, zLoad, lengthToLambda(d1))
    zStub1 = zIn(z0, Inf, lengthToLambda(l1))
    zLoad1 = parallel(zIn1, zStub1)

    zIn2 = zIn(z0, zLoad1, lengthToLambda(d2))
    zStub2 = zIn(z0, Inf, lengthToLambda(l2))
    zLoad2 = parallel(zIn2, zStub2)

    zIn3 = zIn(z0, zLoad2, lengthToLambda(d3))
    zStub3 = zIn(z0, Inf, lengthToLambda(l3))
    zInTotal = parallel(zIn3, zStub3)

    return (zInTotal - 50) / (zInTotal + 50)
end

function stubMinimizingFunction(lengthsVector)
    normf = 0.01 : 0.01 : 2
    Γ = 20log10.(abs.(stubCircuitReflection.(Ref(lengthsVector), normf)))
    return mean(Γ)
end

function stubPlot(lengthsVector)
    normf = 0.01 : 0.01 : 2
    Γ = 20log10.(abs.(stubCircuitReflection.(Ref(lengthsVector), normf)))
    p = Plots.plot(normf, Γ, gridlinewidth=2,
            xlabel="Normalized Frequency" * L"\frac{f}{f_0}", 
            ylabel="Γ in dB")
    display(p)
end

gr()
f0 = 1e9
lambda = 300000000 / f0
upper = lambda * ones(6)
lower = 0.001*lambda .* ones(6) + [0, 0, 0.05*lambda, 0, 0.05*lambda, 0]
constraints = BoxConstraints(lower, upper)
ga = GA(populationSize = 500, selection = tournament(50), crossover=SPX, mutation = PLM())

# Γ = stubMinimizingFunction(rand(Uniform(0.05, 1), 6))
options = Evolutionary.Options(show_trace = true, iterations = 100)

result = Evolutionary.optimize(stubMinimizingFunction, constraints, ga, options)
println("Optimal lengths: ", Evolutionary.minimizer(result))
println("Minimum mean Γ: ", Evolutionary.minimum(result))
stubPlot(Evolutionary.minimizer(result))

summary(result)