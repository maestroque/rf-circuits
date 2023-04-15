
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

function circuit1(f)
    loadLine = 0.2 * f / f0
    shortLine = 0.13 * f / f0
    connectingLine = 0.1 * f / f0


end

N = 201
f0 = 1e9

# Frequencies to be simulated
f = 0 : (4 * f0 / N) : (4 * f0)

z1 = zIn(50, 100 - 79.57im, 0.2) / 50
zShort = zIn(50, 0, 0.13) / 50

z2 = 1 / (1 / z1 + 1 / zShort)

zIn2 = zIn(50, z2, 0.1) / 50

println(zIn2)

zLoadTotal = 1 / zIn2 + 2π * 2.7 * 10e-3 * 50

println(zLoadTotal)

Γ = (zLoadTotal - 50) / (zLoadTotal + 50)

println("Reflection Coefficient Magnitude: ", abs(Γ))