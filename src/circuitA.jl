using Revise
using Plots
using Statistics
using Evolutionary
using Distributions
# includet("transmissionLines.jl")

"""
z0: Characteristic Impedance in Ω
zLoad: Load's Impedance Ω
length: Transmission Line Length measured in fraction of λ

Returns the input impedance of the transmission line
"""
function zIn(z0, zLoad, length)
    if zLoad == Inf
        return -im * z0 * cot(length)
    else
        return z0 * (zLoad + im * z0 * tan(2π * length)) / (z0 + im * zLoad * tan(2π * length))
    end
end

"""
z1, z2: Impedance values in parallel

Returns the equivalent parallel impedance value
"""
function parallel(z1, z2)
    return 1 / (1 / z1 + 1 / z2)
end    

"""
f: Input signal frequency in Hz

Returns the Γ reflection coefficient for the circuit 1.1
for a specific input frequency
"""
function circuitReflection(f)
    z0 = 50
    f0 = 1e9
    zCLoad = 1 / (2π * f * 2 * 10e-12 * im)
    zCIn = 1 / (2π * f * 2.7 * 10e-12 * im)
    loadLine = 0.2 * f / f0
    shortLine = 0.13 * f / f0
    connectingLine = 0.1 * f / f0

    z1 = zIn(z0, 100 + zCLoad, loadLine) 
    zShort = zIn(z0, 0, shortLine) 
    z2 = parallel(z1, zShort)
    zIn2 = zIn(z0, z2, connectingLine) 
    zLoadTotal = parallel(zIn2, zCIn) 
    return (zLoadTotal - 50) / (zLoadTotal + 50)
end

function filterReflection(f)
    length = 0.125 * f / f0
    z0 = 50
    z1 = 98.45
    z2 = 43.6
    z3 = 101.6

    zStub1 = zIn(z1, Inf, length)
    zLoad1 = parallel(zStub1, z0)

    zIn2 = zIn(z3, zLoad1, length)
    zStub2 = zIn(z2, Inf, length)
    zLoad2 = parallel(zIn2, zStub2)

    zIn3 = zIn(z3, zLoad2, length)
    zStub3 = zIn(z1, Inf, length)
    zInTotal = parallel(zIn3, zStub3)

    return (zInTotal - 50) / (zInTotal + 50)

end

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
    plot(normf, Γ)
end

N = 201
f0 = 1e9
z0 = 50

# Frequencies to be simulated
f = 0 : (4 * f0 / N) : (3 * f0)

# Reflection Coefficients for the frequency sweep
# Γ = filterReflection.(f)
# SWR = (1 .+ abs.(Γ)) ./ (1 .- abs.(Γ))

# display(plot(f, 20log10.(abs.(Γ))))
# display(plot(f, SWR))

upper = ones(6)
lower = 0.05 .* ones(6)
constraints = BoxConstraints(lower, upper)
ga = GA(populationSize = 2000)
normf = 0.5 : 0.01 : 1.5

# Γ = stubMinimizingFunction(rand(Uniform(0.05, 1), 6))
result = Evolutionary.optimize(stubMinimizingFunction, constraints, ga)

println("Optimal lengths: ", Evolutionary.minimizer(result))
println("Minimum mean Γ: ", Evolutionary.minimum(result))
stubPlot(Evolutionary.minimizer(result))


