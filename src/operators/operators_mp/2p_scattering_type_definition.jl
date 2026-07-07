# abstract supertype for two-particle scattering operators
abstract type AbstractMP2PScatteringOperator{MPB} <: AbstractMP2POperator{MPB} end

export AbstractMP2PScatteringOperator

# define a MP operator for 2 particle scattering interactions of ELECTRONS
mutable struct MPElectron2PScatteringOperator{
        MPB <: MPBasis{N,SPBS} where {N,SPBS}
    } <: AbstractMP2PScatteringOperator{MPB}

    # the MP basis
    basis :: MPB
    # the matrix representation
    matrix_rep :: SparseMatrixCSC{Complex{Float64}}
    # the list of interacting orbitals ((a,b),(c,d))
    # meaning: b->a, d->c scattering
    interacting_orbitals :: Vector{Tuple{Tuple{Int64,Int64},Tuple{Int64,Int64}}}
    # the prefactor
    prefactor :: Float64

    # custom constructor
    function MPElectron2PScatteringOperator(basis :: MPB, interacting_orbitals :: Vector{Tuple{Tuple{Int64,Int64},Tuple{Int64,Int64}}}, prefactor :: Real) where {
                N, SPBS <: AbstractSPBasisState,
                MPB <: MPBasis{N,SPBS}
            }
        # create a new operator
        op = new{MPB}(basis, spzeros(length(basis), length(basis)), interacting_orbitals, prefactor)
        # recalculate the matrix representation
        recalculate!(op)
        # return the operator
        return op
    end
end
export MPElectron2PScatteringOperator

# define a MP operator for 2 particle scattering interactions of HOLES
mutable struct MPHole2PScatteringOperator{
        MPB <: MPBasis{N,SPBS} where {N,SPBS}
    } <: AbstractMP2PScatteringOperator{MPB}

    # the MP basis
    basis :: MPB
    # the matrix representation
    matrix_rep :: SparseMatrixCSC{Complex{Float64}}
    # the list of interacting orbitals ((a,b),(c,d))
    # meaning: b->a, d->c scattering
    interacting_orbitals :: Vector{Tuple{Tuple{Int64,Int64},Tuple{Int64,Int64}}}
    # the prefactor
    prefactor :: Float64

    # custom constructor
    function MPHole2PScatteringOperator(basis :: MPB, interacting_orbitals :: Vector{Tuple{Tuple{Int64,Int64},Tuple{Int64,Int64}}}, prefactor :: Real) where {
                N, SPBS <: AbstractSPBasisState,
                MPB <: MPBasis{N,SPBS}
            }
        # create a new operator
        op = new{MPB}(basis, spzeros(length(basis), length(basis)), interacting_orbitals, prefactor)
        # recalculate the matrix representation
        recalculate!(op)
        # return the operator
        return op
    end
end
export MPHole2PScatteringOperator


import Base.show
function Base.show(io::IO, op::MPHole2PScatteringOperator{MPB}) where {
            N,
            SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS}
        }
    if haskey(io, :compact)
        print(io, "2-particle (hole-hole) scattering with "*string(length(op.interacting_orbitals))*" (c*a)*(c*a) terms")
    else
        print(io, "2-particle (hole-hole) scattering\n")
        print(io, "Interactions include "*string(length(op.interacting_orbitals))*" (c*a)*(c*a) terms:\n")
        for orbital_pair in op.interacting_orbitals
            orb_1_1 = summary(basis(op).single_particle_basis[orbital_pair[1][1]], "{}")
            orb_1_2 = summary(basis(op).single_particle_basis[orbital_pair[1][2]], "{}")
            orb_2_1 = summary(basis(op).single_particle_basis[orbital_pair[2][1]], "{}")
            orb_2_2 = summary(basis(op).single_particle_basis[orbital_pair[2][2]], "{}")
            print(io, " + (c_"*orb_1_1*"*a_"*orb_1_2*")*(c_"*orb_2_1*"*a_"*orb_2_2*")\n")
        end
        print(io, "Overall prefactor is "*string(op.prefactor)*"\n")
        print(io, "Multi-particle basis contains "*string(length(basis(op)))*" states in total, with "*string(N)*" particles per state\n")
    end
end
function Base.show(io::IO, op::MPElectron2PScatteringOperator{MPB}) where {
            N,
            SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS}
        }
    if haskey(io, :compact)
        print(io, "2-particle (electron-electron) scattering with "*string(length(op.interacting_orbitals))*" (c*a)*(c*a) terms")
    else
        print(io, "2-particle (electron-electron) scattering\n")
        print(io, "Interactions include "*string(length(op.interacting_orbitals))*" (c*a)*(c*a) terms:\n")
        for orbital_pair in op.interacting_orbitals
            orb_1_1 = summary(basis(op).single_particle_basis[orbital_pair[1][1]], "{}")
            orb_1_2 = summary(basis(op).single_particle_basis[orbital_pair[1][2]], "{}")
            orb_2_1 = summary(basis(op).single_particle_basis[orbital_pair[2][1]], "{}")
            orb_2_2 = summary(basis(op).single_particle_basis[orbital_pair[2][2]], "{}")
            print(io, " + (c_"*orb_1_1*"*a_"*orb_1_2*")*(c_"*orb_2_1*"*a_"*orb_2_2*")\n")
        end
        print(io, "Overall prefactor is "*string(op.prefactor)*"\n")
        print(io, "Multi-particle basis contains "*string(length(basis(op)))*" states in total, with "*string(N)*" particles per state\n")
    end
end



##############################################################
#   Interface functions
##############################################################

# obtain the current basis (ELECTRON & HOLE)
function basis(operator :: MPDDOP) :: MPB where {
            N, SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS},
            MPDDOP <: AbstractMP2PScatteringOperator{MPB}
        }
    return operator.basis
end

# obtain the matrix representation (ELECTRON & HOLE)
function matrix_representation(operator :: MPDDOP) :: SparseMatrixCSC{Complex{Float64}}  where {
            N, SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS},
            MPDDOP <: AbstractMP2PScatteringOperator{MPB}
        }
    return operator.matrix_rep .* operator.prefactor
end

# possibly recalculate the matrix representation (ELECTRON & HOLE) (Fallback for non XYZ)
function recalculate!(operator :: MPDDOP, basis_change::Bool=true) where {
            N, SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS},
            MPDDOP <: AbstractMP2PScatteringOperator{MPB}
        }
    @error "currently only recalculation of density-density operator implemented for XYZ basis states, not basis states of type $(SPBS)" stacktrace()
end

# possibly recalculate the matrix representation (ELECTRON & HOLE)
function recalculate!(operator :: MPDDOP, basis_change::Bool=true) where {
            N, SPBS <: Union{SPMSBasisState{BasisStateXYZ}, BasisStateXYZ},
            MPB <: MPBasis{N,SPBS},
            MPDDOP <: AbstractMP2PScatteringOperator{MPB}
        }
    # create new matrix
    operator.matrix_rep = spzeros(Complex{Float64}, length(basis(operator)), length(basis(operator)))
    # iterate over all interacting orbitals
    for orb in operator.interacting_orbitals
        # let it fill its contribution
        fill_orbital_contribution!(operator, orb[1], orb[2])
    end
end



# contribution of one orbital (ELECTRON)
# we want to compute <alpha|c†_a c_b c†_c c_d|beta>, with scattering_1=(a,b), scattering_2=(c,d)
function fill_orbital_contribution!(operator :: MPElectron2PScatteringOperator{MPB}, scattering_1::Tuple{Int64,Int64}, scattering_2::Tuple{Int64,Int64}) where {
            N, SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS}
        }
    # unpack indices
    a, b = scattering_1  # c†_a c_b
    c, d = scattering_2  # c†_c c_d
    # reusable occupation buffer
    occ_buffer = Vector{Int64}(undef, N)
    # loop over beta states containing d 
    for beta  in basis(operator).lookup_sp_states[d]
        state_2 = basis(operator)[beta] # |beta>
        # copy occupation into buffer
        @inbounds for i in 1:N
            occ_buffer[i] = state_2.occupation[i]
        end
        # find position of b -- must be present for c_b to give nonzero (d is guaranteed by lookup_sp_states[d], b must be checked explicitly)
        b_pos = findfirst(==(b), occ_buffer)
        if b_pos === nothing
            continue  # b not in this state -> zero
        end
        # apply c†_c c_d|beta>: replace d -> c
        d_pos = findfirst(==(d), occ_buffer)
        @inbounds occ_buffer[d_pos] = c
        # Pauli: c†_c gives zero if c was already present elsewhere
        if count(==(c), occ_buffer) > 1
            continue
        end
        # apply c†_a c_b: replace b -> a
        @inbounds occ_buffer[b_pos] = a
        # Pauli: c†_a gives zero if a was already present elsewhere
        if count(==(a), occ_buffer) > 1
            continue
        end
        # sort to canonical form, accumulating the full fermionic sign
        s = permutation_sign!(occ_buffer)
        # s==0 means duplicate entries 
        if s == 0
            continue
        end
        # direct Dict lookup to find <alpha|-- occ_buffer is now in canonical
        # (sorted) form, matching the keys stored by resetLookupIndex!
        find_alpha = get(basis(operator).lookup_index, occ_buffer, nothing)
        if find_alpha === nothing
            continue
        end
        alpha = find_alpha[2]
        # fill matrix entry
        operator.matrix_rep[alpha, beta] += s
    end

end

# contribution of one orbital (HOLE)
# we want to compute <alpha|c†_b c_a c†_d c_c |beta> (indices swapped relative to electron case), with scattering_1=(a,b), scattering_2=(c,d)
function fill_orbital_contribution!(operator :: MPHole2PScatteringOperator{MPB}, scattering_1::Tuple{Int64,Int64}, scattering_2::Tuple{Int64,Int64}) where {
            N, SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS}
        }
    # unpack indices -- hole operator is c†_b c_a c†_d c_c
    a, b = scattering_1  # c†_b c_a  
    c, d = scattering_2  # c†_d c_c  
    # reusable occupation buffer
    occ_buffer = Vector{Int64}(undef, N)
    # loop over beta states containing a
    for beta  in basis(operator).lookup_sp_states[a]
        state_2 = basis(operator)[beta] #|beta>
        # copy occupation into buffer
        @inbounds for i in 1:N
            occ_buffer[i] = state_2.occupation[i]
        end
        # find position of c -- must be present for c_c to give nonzero
        c_pos = findfirst(==(c), occ_buffer)
        if c_pos === nothing
            continue  # c not in this state -> zero, c_a c_c|beta>=0
        end
        # apply c†_b c_a|beta>: replace a -> b
        a_pos = findfirst(==(a), occ_buffer)
        @inbounds occ_buffer[a_pos] = b
        # Pauli: c†_b gives zero if b was already present elsewhere
        if count(==(b), occ_buffer) > 1
            continue
        end
        # apply c†_d c_c: replace c -> d
        @inbounds occ_buffer[c_pos] = d
        # Pauli: c†_d gives zero if d was already present elsewhere
        if count(==(d), occ_buffer) > 1
            continue
        end
        # sort to canonical form, accumulating the full fermionic sign
        s = permutation_sign!(occ_buffer)
        # s==0 means duplicate entries
        if s == 0
            continue
        end
        # direct Dict lookup to find <alpha|-- occ_buffer is now in canonical
        # (sorted) form, matching the keys stored by resetLookupIndex!
        find_alpha = get(basis(operator).lookup_index, occ_buffer, nothing)
        if find_alpha === nothing
            continue
        end
        alpha = find_alpha[2]
        # fill matrix entry 
        operator.matrix_rep[alpha, beta] += s
    end
end


export fill_orbital_contribution!
