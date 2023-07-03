using NLsolve
using LaTeXStrings

struct WaveguideMaterial
    dA::Vector{<:Real}
    dB::Vector{<:Real}
    w
end

# Dimensions of a WR-90 waveguide
a = 22.86
b = 10.16

f = 1e9
c = 3e8
η = 120π

# TE10 mode  
fc = c / (2 * a)
# λg = c / (f * sqrt(1 - (fc / f)^2))
Zo = η / sqrt(1 - (fc / f)^2) 

short = [48.4, 73.0, 97.5, 122.2]

material1 = WaveguideMaterial([46.9, 71.5, 96.0, 120.7], [55.6, 80.2, 104.7, 129.4], 1.5)
material2 = WaveguideMaterial([46.7, 71.3, 95.8, 120.5], [52.6, 77.0, 101.7, 126.3], 1.517)

# Shorted waveguide
λg = 2 * sum([short[i + 1] - short[i] for i in range(1, 3)]) / length(short) * 0.001

# Material 1
dAmin1 = 0.001 * sum(material1.dA[i + 1] - short[i] for i in range(1, 3)) / (length(material1.dA))
Zsc1 = -im*Zo*tan(2π / λg * dAmin1)

dBmin1 = 0.001 * sum(material1.dB[i + 1] - short[i] for i in range(1, 3)) / (length(material1.dA))
Zoc1 = -im*Zo*tan(2π / λg * dBmin1)

# Dielectric constant using both measurements
ϵr1_AB = real((fc / f)^2 + η^2 / (Zoc1*Zsc1))

# Using only measurement A 

# Transcedental equation to be solved for x
function g(x)
    2π * f * material1.w * η / c * tan(x[1]) / x[1] - imag(Zsc1)
end

sol = nlsolve(g, [15.0])
x1 = sol.zero[1]

ϵr1_A = (c * x1 / (2\pi * f * material1.w))^2 + (fc / f)^2

# Material 2
dAmin2 = 0.001 * sum(material2.dA[i + 1] - short[i] for i in range(1, 3)) / (length(material2.dA))
Zsc2 = -im*Zo*tan(2π / λg * dAmin2)

dBmin2 = 0.001 * sum(material2.dB[i + 1] - short[i] for i in range(1, 3)) / (length(material2.dA))
Zoc2 = -im*Zo*tan(2π / λg * dBmin2)

# Dielectric constant using both measurements
ϵr2_AB = real((fc / f)^2 + η^2 / (Zoc2*Zsc2))

# Using only measurement A 

# Transcedental equation to be solved for x
function h(x)
    2π * f * material2.w * η / c * tan(x[1]) / x[1] - imag(Zsc2)
end

sol = nlsolve(h, [23.0])
x2 = sol.zero[1]

ϵr2_A = (c * x2 / (2\pi * f * material2.w))^2 + (fc / f)^2

println("Material 1: ϵr = ", ϵr1_A, " (Measurement A) | ", ϵr1_AB, " (Both measurements) ", x1)
println("Material 2: ϵr = ", ϵr2_A, " (Measurement A) | ", ϵr2_AB, " (Both measurements) ", x2)

Plots.plot(0:0.1:30, h.(0:0.1:30), 
            title="Plot of "*L"\frac{2π f w η_0}{c} \frac{tan(x)}{x} - Z_{SC}")