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
    # dimension of the many-particle basis
    D = length(basis(operator))
    # collect (row, col, value) triplets (instead of writing into the sparse matrix one entry at a time)
    Is = Int[]
    Js = Int[]
    exp_vals = Complex{Float64}[]
    # allocate a buffer state
    state_buffer = deepcopy(basis(operator)[1])
    state_buffer.basis_index = -1
    state_buffer.basis_sign  = 0
    # calculate the own matrix elements
    for a in 1:length(operator.basis.single_particle_basis)
    for b in 1:length(operator.basis.single_particle_basis)
        # get the element of the single particle hamiltonian
        op_sp_ab = matrix_rep_sp[a,b]
        # check if nonzero, else continue to next iteration
        if abs(op_sp_ab) <= 1e-8
            continue
        end
        # generate all element contributions to the many body hamiltonian
        for alpha in basis(operator).lookup_sp_states[a]
        for beta  in basis(operator).lookup_sp_states[b]
            # add the expectation with ab to the matrix
            exp = expectation_value_ca!(basis(operator), basis(operator)[alpha], a,b, basis(operator)[beta], state_buffer) * op_sp_ab
            # skip exact zeros so we don't store them as structural nonzeros
            if exp != 0
                push!(Is, alpha)
                push!(Js, beta)
                push!(exp_vals, exp)
            end
        end
        end
    end
    end
    # build the sparse matrix in one pass; combine=+ sums any duplicate
    # (alpha,beta) entries that arose from different (a,b) contributions
    matrix_rep = sparse(Is, Js, exp_vals, D, D, +)
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
