using Revise
using Plots
using PyPlot
using PyCall

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
function antennaArrayFieldIntensity(f, I, d, l, θ, ϕ, N, alternatingCurrents=false)
    λ = 3e8 / f
    k = 2π / λ
    r = 10
    arrayFactor = 1
    for i in range(1, N-1)
        if alternatingCurrents
            arrayFactor += (-1)^(i) * exp(im*k*d*(i)*cos(ϕ)*sin(θ))
        else
            arrayFactor += exp(im*k*d*(i+1)*cos(ϕ)*sin(θ))
        end 
    end

    return abs(im*60*I*exp(-im*k*r)*(cos(k/2 * l * cos(θ)) - cos(k/2 * l)) / (r * sin(θ)) * arrayFactor)
end

function color_function(x, y, z)
    radius = sqrt(x^2 + y^2 + z^2)
    return radius
end

maxSteps = 100
θ = range(0, stop=2π, length=maxSteps)
ϕ = range(0, stop=π, length=maxSteps)
f = 1e9
λ = 3e8 / f
N = 8
d = λ/8
l = λ/2
alternating = true

E = antennaArrayFieldIntensity.(1e9, 1, d, l, 0, ϕ, N)
gr()

display(Plots.plot(ϕ, antennaArrayFieldIntensity.(1e9, 1, d, l, π/2, ϕ, N, alternating), proj=:polar, gridlinewidth=2, title="θ=0"))
display(Plots.plot(θ, antennaArrayFieldIntensity.(1e9, 1, d, l, θ, 0, N, alternating), proj=:polar, gridlinewidth=2, title="φ=0"))

X(r,theta,phi) = r * sin(theta) * sin(phi)
Y(r,theta,phi) = r * sin(theta) * cos(phi)
Z(r,theta,phi) = r * cos(theta)

pyplot()

xs = [X(antennaArrayFieldIntensity(1e9, 1, d, l, theta, phi, N, alternating), theta, phi) for theta in θ, phi in ϕ]
ys = [Y(antennaArrayFieldIntensity(1e9, 1, d, l, theta, phi, N, alternating), theta, phi) for theta in θ, phi in ϕ]
zs = [Z(antennaArrayFieldIntensity(1e9, 1, d, l, theta, phi, N, alternating), theta, phi) for theta in θ, phi in ϕ]

# Initiate the radial color map
R = color_function.(xs, ys, zs)
R[1, :] = 0*ones(size(zs[1,:]))
R = R ./ maximum(R)

# Plot the 3D radiation pattern
surface(xs, ys, zs, fill_z = R)
