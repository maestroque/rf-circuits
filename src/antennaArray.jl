using Revise
using Plots
using PyPlot
using PyCall

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

"""
## Arguements
- `f`: Frequency
- `I`: Current 
- `d`: Distance between dipoles
- `θ`: 
- `ϕ`:
- `N`: Number of dipoles
- `alternatingCurrents`: Whether all currents are the same polarity or alternating
"""
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

function color_function(x, y, z)
    radius = sqrt(x^2 + y^2 + z^2)
    return radius
end

maxSteps = 100
f = 1e9
I = 1
θ = range(0, stop=2π, length=maxSteps)
ϕ = range(0, stop=π, length=maxSteps)
a = AntennaArray(f, I, 0.5, 0.5)


# E = antennaArrayFieldIntensity.(Ref(a), 0, ϕ)
gr()

display(Plots.plot(ϕ, antennaArrayFieldIntensity.(Ref(a), π/2, ϕ), proj=:polar, gridlinewidth=2, title="θ=0"))
display(Plots.plot(θ, antennaArrayFieldIntensity.(Ref(a), θ, 0), proj=:polar, gridlinewidth=2, title="φ=0"))

X(r,theta,phi) = r * sin(theta) * sin(phi)
Y(r,theta,phi) = r * sin(theta) * cos(phi)
Z(r,theta,phi) = r * cos(theta)

pyplot()

xs = [X(antennaArrayFieldIntensity(a, theta, phi), theta, phi) for theta in θ, phi in ϕ]
ys = [Y(antennaArrayFieldIntensity(a, theta, phi), theta, phi) for theta in θ, phi in ϕ]
zs = [Z(antennaArrayFieldIntensity(a, theta, phi), theta, phi) for theta in θ, phi in ϕ]

# Initiate the radial color map
R = color_function.(xs, ys, zs)
R[1, :] = 0*ones(size(zs[1,:]))
R = R ./ maximum(R)

# Plot the 3D radiation pattern
surface(xs, ys, zs, fill_z = R)
