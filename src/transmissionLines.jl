"""
## Arguements
- `z0`: Characteristic Impedance in Ω
- `zLoad`: Load's Impedance Ω
- `length`: Transmission Line Length measured in fraction of λ

Returns the input impedance of the transmission line
"""
function zIn(z0::Number, zLoad::Number, length::Real)
    if !isfinite(zLoad)
        return -im * z0 * cot(2π * length)
    else
        return z0 * (zLoad + im * z0 * tan(2π * length)) / 
        (z0 + im * zLoad * tan(2π * length))
    end
end

"""
## Arguements
- `z1::Complex`, `z2::Complex`: Impedance values

Returns the equivalent parallel impedance value
"""
function parallel(z1::Complex, z2::Complex)
    return 1 / (1 / z1 + 1 / z2)
end    