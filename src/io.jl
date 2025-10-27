import HDF5

"""
	pellet2hdf(pelt::Pellet, filename::AbstractString)

Write a Pellet struct to an HDF5 file. Excludes private fields (e.g., `_pool`).

# Arguments
- `pelt::Pellet`: Pellet struct to write
- `filename::AbstractString`: Output HDF5 filename (must end with .h5)

# Notes
- Uses IMAS.imas2hdf to serialize the `properties` field
- Automatically converts Symbol fields to Strings for HDF5 compatibility
"""
function pellet2hdf(pelt::Pellet, filename::AbstractString)
	@assert endswith(filename, ".h5") "filename must end with .h5"


	HDF5.h5open(filename, "w") do fid
		# Write properties using IMAS converter
		gparent = HDF5.create_group(fid, "properties")
		IMAS.imas2hdf(pelt.properties, gparent)

		# Write all other fields (excluding properties and private_fields)
		for field in fieldnames(typeof(pelt))
			if field in (:properties, private_fields...)
				continue
			end

			value = getfield(pelt, field)

			# Convert Symbol to String for HDF5 compatibility
			if value isa Symbol
				fid[string(field)] = string(value)
			else
				fid[string(field)] = value
			end
		end
	end
end



"""
	hdf2pellet(filename::AbstractString) -> Pellet

Read a Pellet struct from an HDF5 file created by `pellet2hdf`.

# Arguments
- `filename::AbstractString`: Input HDF5 filename (must end with .h5)

# Returns
- `Pellet{FT}`: Reconstructed Pellet struct with inferred floating point type

# Notes
- Automatically infers floating point type (Float32/Float64) from stored data
- Uses IMAS.hdf2imas to deserialize the `properties` field
- Automatically converts String fields back to Symbols where appropriate
- Converts arrays to StaticArrays.MVector for fixed-size vector fields
"""
function hdf2pellet(filename::AbstractString)

	fid = HDF5.h5open(filename, "r") 

	# Determine the floating point type from stored data
	FT = typeof(fid["t"][]) 

	# Create Pellet instance
	pelt = Pellet{FT}()

	# Read properties using IMAS converter
	pelt.properties = IMAS.pellets__time_slice___pellet()
	IMAS.hdf2imas(fid["properties"], pelt.properties; 
				show_warnings=false, skip_non_coordinates=false, error_on_missing_coordinates=false)

	# Read all other fields (excluding properties and private_fields)
	for field in fieldnames(typeof(pelt))
		if field in (:properties, private_fields...)
			continue
		end

		value = fid[string(field)][]

		# Convert String back to Symbol if necessary
		if fieldtype(Pellet, field) === Symbol
			setfield!(pelt, field, Symbol(value))
		elseif fieldtype(Pellet, field) <: StaticArrays.MVector
			setfield!(pelt, field, StaticArrays.MVector{length(value), FT}(value))
		else
			setfield!(pelt, field, value)
		end
	end
	
	close(fid)

	return pelt
end
