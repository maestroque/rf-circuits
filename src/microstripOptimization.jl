function stubCircuitReflection(lengthsVector, f)
    d1 = lengthsVector[1]
    l1 = lengthsVector[2]
    d2 = lengthsVector[3]
    l2 = lengthsVector[4]
    d3 = lengthsVector[5]
    l3 = lengthsVector[6]

    zLoad = 120 - 80im
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
    normf = 0.5 : 0.01 : 1.5
    Γ = 20log10.(abs.(stubCircuitReflection.(Ref(lengthsVector), normf)))
    return mean(Γ)
end

function stubPlot(lengthsVector)
    normf = 0.5 : 0.01 : 1.5
    Γ = 20log10.(abs.(stubCircuitReflection.(Ref(lengthsVector), normf)))
    display(plot(normf, Γ))
end

# upper = ones(6)
# lower = 0.001 .* ones(6) + [0, 0, 0.05, 0, 0.05, 0]
# constraints = BoxConstraints(lower, upper)
# ga = GA(populationSize = 200, selection = tournament(20), crossover=SPX, mutation = PLM())

# # Γ = stubMinimizingFunction(rand(Uniform(0.05, 1), 6))
# options = Evolutionary.Options(show_trace = true, iterations = 100)

# result = Evolutionary.optimize(stubMinimizingFunction, constraints, ga, options)
# println("Optimal lengths: ", Evolutionary.minimizer(result))
# println("Minimum mean Γ: ", Evolutionary.minimum(result))
# stubPlot(Evolutionary.minimizer(result))

# summary(result)