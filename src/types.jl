
import StaticArrays

"""
Preallocated buffer pool to reduce allocations
"""
mutable struct AllocationPool{T<:Real}
    nsource_buffer::Vector{T}
    nsource_max_len::Int
end

function AllocationPool(::Type{T}, surfaces::Vector{<:IMAS.AbstractFluxSurface}) where {T<:Real}
    max_len = maximum(length(surf.r) for surf in surfaces)
    return AllocationPool{T}(Vector{T}(undef, max_len), max_len)
end

mutable struct Pellet{FT<:AbstractFloat}
    properties::IMAS.pellets__time_slice___pellet
    Btdep::Bool
    update_plasma::Bool
    drift_model::Symbol
    time::Vector{FT}
    t::FT  # initial time
    Bt::Vector{FT}

    velocity_vector::StaticArrays.MVector{3, FT}

    r::Vector{FT}
    R_drift::Vector{FT}
    z::Vector{FT}
    x::Vector{FT}
    y::Vector{FT}
    ρ::Vector{FT}
    Te::Vector{FT}
    ne::Vector{FT}
    radius::Vector{FT}
    ablation_rate::Vector{FT}
    density_source::Matrix{FT}
    temp_drop::Vector{FT}

    # private field for pre-allocation
    _pool::AllocationPool{FT}

    # Empty constructor for reading from file
    function Pellet{FT}() where {FT<:AbstractFloat}
        new{FT}()
    end

    # Full constructor for outer constructor compatibility
    function Pellet{FT}(properties, Btdep, update_plasma, drift_model, time, t, Bt,
                        velocity_vector, r, R_drift, z, x, y, ρ, Te, ne, radius,
                        ablation_rate, density_source, temp_drop, _pool) where {FT<:AbstractFloat}
        new{FT}(properties, Btdep, update_plasma, drift_model, time, t, Bt,
                velocity_vector, r, R_drift, z, x, y, ρ, Te, ne, radius,
                ablation_rate, density_source, temp_drop, _pool)
    end
end

const private_fields = (:_pool, )
