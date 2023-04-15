using Plots

#=
z0: Characteristic Impedance in Ω
zLoad: Load's Impedance Ω
length: Transmission Line Length measured in fraction of λ
=#
function zIn(z0, zLoad, length)
    if length == 0.25
        return z0 ^ 2 / zLoad
    else
        return z0 * (zLoad + im * z0 * tan(2π * length)) / (z0 + im * zLoad * tan(2π * length))
    end
end

function circuitReflection(f)
    z0 = 50
    loadLine = 0.2 * f / f0
    shortLine = 0.13 * f / f0
    connectingLine = 0.1 * f / f0

    z1 = zIn(z0, 100 - 79.57im, loadLine) / z0
    zShort = zIn(z0, 0, shortLine) / z0
    z2 = 1 / (1 / z1 + 1 / zShort)
    zIn2 = zIn(z0, z2 * z0, connectingLine) / z0
    zLoadTotal = 1 / zIn2 + 2π * 2.7 * 10e-3 * 50
    return (zLoadTotal - 50) / (zLoadTotal + 50)
end

N = 201
f0 = 1e9
z0 = 50

# Frequencies to be simulated
f = 0 : (4 * f0 / N) : (4 * f0)

# Reflection Coefficients for the frequency sweep
Γ = circuitReflection.(f)

plot(f, 20log10.(abs.(Γ)))
