################################################################################
#   Type definition
################################################################################


# Type Definition of SumOperator
"""
    mutable struct SumOperator{B, O1<:AbstractOperator{B}, O2<:AbstractOperator{B}} <: AbstractOperator{B}

This object defines the sum operator.

# Fields

- `basis :: B`, the basis;
- `op_1 :: O1`, `op_2 :: O2`, the two contained operators.

"""
mutable struct SumOperator{B, O1<:AbstractOperator{B}, O2<:AbstractOperator{B}} <: AbstractOperator{B}
    # the basis
    basis :: B
    # the two contained operators
    op_1 :: O1
    op_2 :: O2

    # Custom constructor (without explicit matrix rep)
    function SumOperator(op_1::O1, op_2::O2) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O1<:AbstractOperator{B}, O2<:AbstractOperator{B}}
        # check that they have the same basis
        @assert basis(op_1) == basis(op_2)
        # construct new operator
        op = new{B, O1, O2}(basis(op_1), op_1, op_2)
        # return the operator
        return op
    end
end

# export operator type
export  SumOperator

import Base.show
function Base.show(io::IO, op::SumOperator{B, O1, O2}) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O1<:AbstractOperator{B}, O2<:AbstractOperator{B}}
    if haskey(io, :compact)
        show(io, op.op_1)
        print(io, "\n+ ")
        show(io, op.op_2)
    else
        print(io, "Operator sum of the following operators:\n  ")
        show(IOContext(io, :compact=>true), op.op_1)
        print(io, "\n+ ")
        show(IOContext(io, :compact=>true), op.op_2)
        print(io, "\n")
        printBasisInformation(io, op)
    end
end


function printBasisInformation(io::IO, op::SumOperator{B, O1, O2}) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O1<:AbstractOperator{B}, O2<:AbstractOperator{B}}
    print(io, "Basis contains "*string(length(basis(op)))*" states in total, states are of type "*string(BS)*"\n")
end

function printBasisInformation(io::IO, op::SumOperator{B, O1, O2}) where {BS<:AbstractSPSSBasisState, B<:AbstractBasis{BS}, O1<:AbstractOperator{B}, O2<:AbstractOperator{B}}
    print(io, "Basis is single-particle / single-site basis and contains "*string(length(basis(op)))*" states of type "*string(BS)*"\n")
end
function printBasisInformation(io::IO, op::SumOperator{B, O1, O2}) where {SPBS<:AbstractSPSSBasisState, BS<:SPMSBasisState{SPBS}, B<:AbstractBasis{BS}, O1<:AbstractOperator{B}, O2<:AbstractOperator{B}}
    print(io, "Basis is single-particle / multi-site basis and contains "*string(length(basis(op)))*" states with SPSS type "*string(SPBS)*"\n")
    print(io, "Basis contains sites "*string(unique([b.site for b in basis(op)]))*"\n")
end
function printBasisInformation(io::IO, op::SumOperator{B, O1, O2}) where {N, MPBS<:MPBasisState{N}, SPBS<:AbstractSPSSBasisState, B<:MPBasis{MPBS,SPBS}, O1<:AbstractOperator{B}, O2<:AbstractOperator{B}}
    print(io, "Basis is multi-particle / single-site basis and contains "*string(length(basis(op)))*" states with "*string(N)*" particles each\n")
    print(io, "SPSS basis contains "*string(length(basis(op.single_particle_basis)))*" states of type "*string(SPBS)*"\n")
end
function printBasisInformation(io::IO, op::SumOperator{B, O1, O2}) where {N, MPBS<:MPBasisState{N}, SPSSBS<:AbstractSPSSBasisState, SPBS<:SPMSBasisState{SPSSBS}, B<:MPBasis{N,SPBS}, O1<:AbstractOperator{B}, O2<:AbstractOperator{B}}
    print(io, "Basis is multi-particle / multi-site basis and contains "*string(length(basis(op)))*" states with "*string(N)*" particles each\n")
    print(io, "SPSS basis contains "*string(length(basis(op).single_particle_basis))*" states of type "*string(SPSSBS)*" on sites "*string(unique([b.site for b in basis(op).single_particle_basis]))*"\n")
end



################################################################################
#   Interface functions
################################################################################

# obtain the current basis
function basis(operator :: SumOperator{B, O1, O2}) :: B where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O1<:AbstractOperator{B}, O2<:AbstractOperator{B}}
    return operator.basis
end

# calculate the matrix representation
function matrix_representation(operator :: SumOperator{B, O1, O2}) :: SparseMatrixCSC{Complex{Float64}} where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O1<:AbstractOperator{B}, O2<:AbstractOperator{B}}
    # return sum of matrix representations
    return matrix_representation(operator.op_1) .+ matrix_representation(operator.op_2)
end

# possibly recalculate the matrix representation
function recalculate!(operator :: SumOperator{B, O1, O2}, basis_change::Bool=true)  where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O1<:AbstractOperator{B}, O2<:AbstractOperator{B}}
    # recalculate contained operators
    recalculate!(operator.op_1, basis_change)
    recalculate!(operator.op_2, basis_change)
end

# set a parameter (returns (found parameter?, changed matrix?))
function set_parameter!(operator :: SumOperator{B, O1, O2}, parameter :: Symbol, value; print_result::Bool=false, recalculate::Bool=true, kwargs...) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O1<:AbstractOperator{B}, O2<:AbstractOperator{B}}
    # check if it can be set in operator 1 or 2
    found_param_1, changed_matrix_1 = set_parameter!(operator.op_1, parameter, value, print_result=print_result, recalculate=recalculate; kwargs...)
    found_param_2, changed_matrix_2 = set_parameter!(operator.op_2, parameter, value, print_result=print_result, recalculate=recalculate; kwargs...)
    return (found_param_1 || found_param_2, changed_matrix_1 || changed_matrix_2)
end

# get a parameter (returns (found parameter?, parameter value or nothing))
function get_parameter(operator :: SumOperator{B, O1, O2}, parameter :: Symbol; print_result::Bool=false, kwargs...) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O1<:AbstractOperator{B}, O2<:AbstractOperator{B}}
    # check if parameter can be obtained
    param_1 = get_parameter(operator.op_1, parameter, print_result=print_result; kwargs...)
    param_2 = get_parameter(operator.op_2, parameter, print_result=print_result; kwargs...)
    # give warning if paramter found in both
    if param_1 != nothing && param_2 != nothing
        if print_result
            println("WARNING, parameter :$(parameter) found in both sub operators!!! Returning value of operator 1")
        end
        return param_1
    elseif param_1 != nothing
        if print_result
            println("Parameter :$(parameter) found in sub operator 1, value $(param_1)")
        end
        return param_1
    elseif param_2 != nothing
        if print_result
            println("Parameter :$(parameter) found in sub operator 2, value $(param_2)")
        end
        return param_2
    else
        if print_result
            println("Parameter :$(parameter) not found in any sub operator")
        end
    end
end

# get a parameter (returns (found parameter?, parameter value or nothing))
function get_parameters(operator :: SumOperator{B, O1, O2}; kwargs...) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O1<:AbstractOperator{B}, O2<:AbstractOperator{B}}
    # check if parameter can be obtained
    params_1 = get_parameters(operator.op_1; kwargs...)
    params_2 = get_parameters(operator.op_2; kwargs...)
    params = Symbol[]
    append!(params, params_1)
    append!(params, params_2)
    return unique(sort(params))
end




################################################################################
#   (+) functions
################################################################################

# import of base function
import Base.(+)

# Define the (+) function for any operators with same basis type
function +(op_1 :: O1, op_2 :: O2) :: SumOperator{B, O1, O2}  where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O1<:AbstractOperator{B}, O2<:AbstractOperator{B}}
    # build a new operator and return it
    return SumOperator(op_1, op_2)
end

# Define the (+) function for any operators with different single particle single site basis types
#function +(op_1 :: O1, op_2 :: O2)  where {BS1<:AbstractBasisState, B1<:AbstractBasis{BS1}, BS2<:AbstractBasisState, B2<:AbstractBasis{BS2}, O1<:AbstractOperator{B1}, O2<:AbstractOperator{B2}}
#    # return projected into basis 1
#    SumOperator(op_1, ProjectorOperator(op_2, basis(op_1)))
#end