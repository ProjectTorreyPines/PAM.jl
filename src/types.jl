
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
const public_fields = filter(x -> x ∉ private_fields, fieldnames(Pellet))


# Equality operators

"""
Programming equality: NaN == NaN, missing == missing. Always returns Bool.
"""
function Base.isequal(p1::Pellet, p2::Pellet)
	for field in public_fields
		isequal(getfield(p1, field), getfield(p2, field)) || return false
	end
	return true
end

"""
Mathematical equality: NaN ≠ NaN, missing propagates. Returns Bool or Missing.
"""
function Base.:(==)(p1::Pellet, p2::Pellet)
	for field in public_fields
		field_eq = getfield(p1, field) == getfield(p2, field)
		ismissing(field_eq) && return missing
		field_eq || return false
	end
	return true
end

"""
Approximate equality with tolerance. Useful for comparing simulation results.
Supports rtol and atol keyword arguments (see Base.isapprox).
Non-numeric fields (Symbol, Bool) are compared with exact equality.
"""
function Base.isapprox(p1::Pellet, p2::Pellet; kwargs...)
	for field in public_fields
		a = getfield(p1, field)
		b = getfield(p2, field)

		# Non-numeric types use exact equality
		if a isa Union{Symbol, Bool}
			isequal(a, b) || return false
		else
			# Numeric types use approximate equality
			isapprox(a, b; kwargs...) || return false
		end
	end
	return true
end

"""
	diff(p1::Pellet, p2::Pellet; verbose=true) -> Bool

Compare two Pellets and print differences. Returns true if different, false if equal.
Set `verbose=false` for fast check without output.
"""
function Base.diff(p1::Pellet, p2::Pellet; verbose::Bool=true)
	has_diff = false
	for field in public_fields
		if !isequal(getfield(p1, field), getfield(p2, field))
			verbose ? println("Field \"$(field)\" differs") : return true
			has_diff = true
		end
	end
	return has_diff
end

