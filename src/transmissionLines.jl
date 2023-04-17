
"""
z0: Characteristic Impedance in Ω
zLoad: Load's Impedance Ω
length: Transmission Line Length measured in fraction of λ

Returns the input impedance of the transmission line
"""
function zIn(z0, zLoad, length)
    if length == 0.25
        return z0 ^ 2 / zLoad
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