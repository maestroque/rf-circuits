using Revise
using Plots
using Statistics
using Evolutionary
using Distributions
includet("transmissionLines.jl")

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

N = 400
f0 = 1e9
z0 = 50

# Frequencies to be simulated
f = 0 : (4 * f0 / N) : (3 * f0)

# Reflection Coefficients for the frequency sweep
Γ = filterReflection.(f)
SWR = (1 .+ abs.(Γ)) ./ (1 .- abs.(Γ))
SWR[SWR .> 10] .= 10

reflectionPlot = plot(
                    f ./ 10e8, 20log10.(abs.(Γ)), gridlinewidth=2,
                    xlabel="Frequency in GHz", ylabel="Γ in dB"
                    )

swrPlot = plot(
                f ./ 10e8, SWR, gridlinewidth=2,
                xlabel="Frequency in GHz", ylabel="SWR")

display(reflectionPlot)
display(swrPlot)