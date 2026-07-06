# define a generalized SP operator type for MP basis states
mutable struct MPGeneralizedSPOperator{
    SPBS <: AbstractSPBasisState,
    MPB <: MPBasis{N,SPBS} where {N},
    SPO <: AbstractSPOperator{SPBasis{SPBS}}
} <: AbstractMP1POperator{MPB}

# the MP basis
basis :: MPB
# the single particle operator
operator :: SPO

# custom constructor
function MPGeneralizedSPOperator(basis :: MPB, operator :: SPO) where {
            N,
            SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS},
            SPO <: AbstractSPOperator{SPBasis{SPBS}}
        }
    # create a new operator
    op = new{SPBS, MPB, SPO}(basis, operator)
    # return the operator
    return op
end
end
export MPGeneralizedSPOperator




import Base.show
function Base.show(io::IO, op::MPGeneralizedSPOperator{SPBS, MPB, SPO}) where {
        N,
        SPBS <: AbstractSPBasisState,
        MPB <: MPBasis{N,SPBS},
        SPO <: AbstractSPOperator{SPBasis{SPBS}}
    }
if haskey(io, :compact)
    print(io, "1-particle ")
    show(io, op.operator)
else
    print(io, "Multi-particle ")
    show(io, op.operator)
    print(io, "Multi-particle basis contains "*string(length(basis(op)))*" states in total, with "*string(N)*" particles per state\n")
end
end




##############################################################
#   Interface functions
##############################################################

# obtain the current basis
function basis(operator :: MPGeneralizedSPOperator{SPBS, MPB, SPO}) :: MPB where {
        N,
        SPBS <: AbstractSPBasisState,
        MPB <: MPBasis{N,SPBS},
        SPO <: AbstractSPOperator{SPBasis{SPBS}}
    }
return operator.basis
end

# calculate the matrix representation
function matrix_representation(operator :: MPGeneralizedSPOperator{SPBS, MPB, SPO}) :: SparseMatrixCSC{Complex{Float64}} where {
        N,
        SPBS <: AbstractSPBasisState,
        MPB <: MPBasis{N,SPBS},
        SPO <: AbstractSPOperator{SPBasis{SPBS}}
    }
    # get matrix representation
    matrix_rep_sp = matrix_representation(operator.operator)
    # get the important matrix elements
    relevant_sp = map(x->abs(x)>1e-8, matrix_rep_sp)
    # create new matrix
    matrix_rep = spzeros(Complex{Float64}, length(basis(operator)), length(basis(operator)))
    # reusable buffer for the occupation vector -- avoids one allocation per inner-loop iteration
    occ_buffer = Vector{Int64}(undef, N)
    # calculate the matrix elements
    for a in 1:length(operator.basis.single_particle_basis)
    for b in 1:length(operator.basis.single_particle_basis)
        # check if relevant
        if !relevant_sp[a,b]
            continue
        end
        # get the element of the single particle hamiltonian
        op_sp_ab = matrix_rep_sp[a,b]
        # generate all element contributions to the many body hamiltonian (in 1 loop): op_sp_ap*<alpha|c_a^\dagger c_b|beta>
        # beta are all the mp states that include sp state b
        for beta in basis(operator).lookup_sp_states[b]
            state_2 = basis(operator)[beta] # |beta>
            # copy occupation into buffer
            # (lookup_sp_states[b] guarantees b is in state_2.occupation)
            @inbounds for i in 1:N
                occ_buffer[i] = state_2.occupation[i]
            end
            # now replace b with a (remove b -> c_b|beta>)
            b_pos = findfirst(==(b), occ_buffer)
            @inbounds occ_buffer[b_pos] = a
            # Pauli exclusion: if a was already present, the result is zero (c_a^\dagger c_b|beta>=0)
            # check cheaply by counting occurrences of a after replacement, if already present skip to next iteration
            if count(==(a), occ_buffer) > 1
                continue
            end
            # sort occ_buffer into canonical form, getting the fermionic sign (c_a^\dagger c_b|beta>=s*|sorted_occ>)
            # permutation_sign! sorts in place and returns the sign in one pass
            s = permutation_sign!(occ_buffer)
            # s==0 means duplicate entries (already caught above, but guard anyway)
            if s == 0
                continue
            end
            # direct Dict lookup to find <alpha|-- O(1) since occ_buffer is now in canonical
            # (sorted) form, matching the keys stored by resetLookupIndex!
            find_alpha = get(basis(operator).lookup_index, occ_buffer, nothing)  # find_alpha = (sign, index); sign from Dict is always 1 for canonical states
            if find_alpha === nothing
                continue
            end
            alpha = find_alpha[2]
            # the full fermionic sign is s from the sort above, matrix element 
            v = s * op_sp_ab # v is always nonzero here (s = ±1, op_sp_ab passed relevant_sp filter)
            # check for floating point cases before saving matrix element
            if v != 0
                matrix_rep[alpha, beta] += v
            end
        end
    end
    end
    #return matrix representation
    return matrix_rep
end


# possibly recalculate the matrix representation
function recalculate!(operator :: MPGeneralizedSPOperator{SPBS, MPB, SPO}, basis_change::Bool=true) where {
        N,
        SPBS <: AbstractSPBasisState,
        MPB <: MPBasis{N,SPBS},
        SPO <: AbstractSPOperator{SPBasis{SPBS}}
    }
    if basis_change
        # reset the basis in the single particle operator
        operator.operator.basis = operator.basis.single_particle_basis
        # let operator recalculate
        recalculate!(operator.operator)
    else
        # let operator recalculate
        recalculate!(operator.operator)
    end
end


# set a parameter (returns (found parameter?, changed matrix?))
function set_parameter!(operator :: MPGeneralizedSPOperator{SPBS, MPB, SPO}, parameter :: Symbol, value; print_result::Bool=true, recalculate::Bool=true, kwargs...) where {
            N,
            SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS},
            SPO <: AbstractSPOperator{SPBasis{SPBS}}
        }
    # pass on to contained operator
    found_param, changed_matrix = set_parameter!(operator.operator, parameter, value, print_result=print_result, recalculate=recalculate; kwargs...)
    return (found_param, changed_matrix)
end

# get a parameter (returns (found parameter?, parameter value or nothing))
function get_parameter(operator :: MPGeneralizedSPOperator{SPBS, MPB, SPO}, parameter :: Symbol; print_result::Bool=true, site::Union{Int64, Symbol}=-1, kwargs...) where {
        N,
        SPBS <: AbstractSPBasisState,
        MPB <: MPBasis{N,SPBS},
        SPO <: AbstractSPOperator{SPBasis{SPBS}}
    }
# pass on to containing operator
return get_parameter(operator.operator, parameter, print_result=print_result; site=site, kwargs...)
end

# get a parameter (returns (found parameter?, parameter value or nothing))
function get_parameters(operator :: MPGeneralizedSPOperator{SPBS, MPB, SPO}; site::Union{Int64, Symbol}=-1, kwargs...) where {
        N,
        SPBS <: AbstractSPBasisState,
        MPB <: MPBasis{N,SPBS},
        SPO <: AbstractSPOperator{SPBasis{SPBS}}
    }
# pass on to containing operator
return get_parameters(operator.operator; site=site, kwargs...)
end
