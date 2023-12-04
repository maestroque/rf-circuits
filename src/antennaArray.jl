using Revise
using Plots
using PlotlyJS

mutable struct AntennaArray
    λ::Float64
    I::Float64
    d::Float64
    l::Float64
    N::Int
    alternatingCurrents::Bool

    """
    ## Arguements
    - `f`: Frequency
    - `I`: Current 
    - `dλ`: Distance between dipoles (fraction of λ)
    - `lλ`: Length of the dipole (fraction of λ)
    - `N`: Number of dipoles
    - `alternatingCurrents`: Whether all currents are the same polarity or alternating
    """
    function AntennaArray(f, I, dλ, lλ, N = 1, alternatingCurrents = false)
        λ = 3e8 / f
        d = dλ * λ
        l = lλ * λ
        new(λ, I, d, l, N, alternatingCurrents)
    end
end

function antennaArrayFieldIntensity(a::AntennaArray, θ, ϕ)
    k = 2π / a.λ
    r = 10
    arrayFactor = 1
    for i in range(1, a.N-1)
        if a.alternatingCurrents
            arrayFactor += (-1)^(i) * exp(im * k * a.d * (i) * cos(ϕ) * sin(θ))
        else
            arrayFactor += exp(im * k * a.d * (i + 1) * cos(ϕ) * sin(θ))
        end 
    end

    return abs(im * 60 * a.I * exp(-im * k * r) * (cos(k / 2 * a.l * cos(θ)) - cos(k / 2 * a.l)) / (r * sin(θ)) * arrayFactor)
end

function directivity(a::AntennaArray)
    maxSteps = 100
    θ = range(2π/maxSteps, stop=2π, length=maxSteps)
    ϕ = range(π/maxSteps, stop=π, length=maxSteps)
    η0 = 120*π

    E = antennaArrayFieldIntensity.(Ref(a), θ, ϕ)
    P = (abs.(E)).^2 ./ (2 * η0)

    # Pav = sum(P .* sin.(θ) .* 2π/50 * π/50)
    Pav = 0
    for p ∈ ϕ
        for t ∈ θ
            El = antennaArrayFieldIntensity(a, p, t)
            Pl = (abs(El))^2 / (2 * η0)
            Pav += Pl * sin(t) * 2π/50 * π/50
            println(Pav)
        end
    end
    println("Total Average", Pav)
    println(maximum(P))
    return 10log10(4π * maximum(P) / Pav)

end

function color_function(x, y, z)
    radius = sqrt(x^2 + y^2 + z^2)
    return radius
end

maxSteps = 100
f = 1e9
I = 1
θ = range(0, stop=2π, length=maxSteps)
ϕ = range(0, stop=π, length=maxSteps)
a = AntennaArray(f, I, 0.25, 0.5, 4, true)

plotlyjs()

display(Plots.plot(ϕ, antennaArrayFieldIntensity.(Ref(a), π/2, ϕ), proj=:polar, gridlinewidth=2, title="θ=0"))
display(Plots.plot(θ, antennaArrayFieldIntensity.(Ref(a), θ, 0), proj=:polar, gridlinewidth=2, title="φ=0"))

X(r,theta,phi) = r * sin(theta) * sin(phi)
Y(r,theta,phi) = r * sin(theta) * cos(phi)
Z(r,theta,phi) = r * cos(theta)

xs = [X(antennaArrayFieldIntensity(a, theta, phi), theta, phi) for theta in θ, phi in ϕ]
ys = [Y(antennaArrayFieldIntensity(a, theta, phi), theta, phi) for theta in θ, phi in ϕ]
zs = [Z(antennaArrayFieldIntensity(a, theta, phi), theta, phi) for theta in θ, phi in ϕ]

# Plot the 3D radiation pattern
PlotlyJS.plot(PlotlyJS.surface(x=xs, y=ys, z=zs, surfacecolor=@. xs^2 + ys^2 + zs^2))