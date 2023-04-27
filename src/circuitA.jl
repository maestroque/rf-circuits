using Revise
using Plots
using Statistics
using Evolutionary
using Distributions
includet("transmissionLines.jl")

"""
f: Input signal frequency in Hz

Returns the Γ reflection coefficient for the circuit 1.1
for a specific input frequency
"""
function circuitReflection(f)
    z0 = 50
    f0 = 1e9
    zRLoad = 100
    zCLoad = 1 / (2π * f * 2 * 10e-12 * im)
    zLoad = zRLoad + zCLoad
    zCIn = 1 / (2π * f * 2.7 * 10e-12 * im)
    loadLine = 0.2 * f / f0
    shortLine = 0.13 * f / f0
    connectingLine = 0.1 * f / f0

    z1 = zIn(z0, zLoad, loadLine) 
    zShort = zIn(z0, 0, shortLine) 
    z2 = parallel(z1, zShort)
    zIn2 = zIn(z0, z2, connectingLine) 
    zLoadTotal = parallel(zIn2, zCIn) 
    return (zLoadTotal - 50) / (zLoadTotal + 50)
end

N = 201
f0 = 1e9
z0 = 50

# Frequencies to be simulated
f = 0 : (4 * f0 / N) : (2 * f0)

# Reflection Coefficients for the frequency sweep
Γ = circuitReflection.(f)
SWR = (1 .+ abs.(Γ)) ./ (1 .- abs.(Γ))

display(plot(f, 20log10.(abs.(Γ))))
display(plot(f, SWR))


println(abs(circuitReflection(0.3e9)))
