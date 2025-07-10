# Type Definition of LS / SpinOrbit Operator

"""
    mutable struct ScalarProductOperator{B, O<:AbstractOperator{B}} <: AbstractOperator{B}

This object defines the scalar product operator.

# Fields

- `basis :: B`, the basis;
- `factor :: Complex{Float64}`, the scalar factor;
- `op :: O`, the contained operator.

"""
mutable struct ScalarProductOperator{B, O<:AbstractOperator{B}} <: AbstractOperator{B}
    # the basis
    basis :: B
    # the scalar factor
    factor :: Complex{Float64}
    # the contained operator
    op :: O

    # Custom constructor (without explicit matrix rep)
    function ScalarProductOperator(factor::Number, op::O) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
        # construct new operator
        operator = new{B, O}(basis(op), factor, op)
        # return the operator
        return operator
    end
end

# export operator type
export ScalarProductOperator

import Base.show
function Base.show(io::IO, op::ScalarProductOperator{B, O}) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    if haskey(io, :compact)
        print(io, "($(op.factor)) * {")
        show(io, op.op)
        print(io, "}")
    else
        print(io, "$(op.factor) times ")
        show(io, op.op)
    end
end




################################################################################
#   Interface functions
################################################################################

# obtain the current basis
function basis(operator :: ScalarProductOperator{B, O}) :: B where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    return operator.basis
end

# obtain the matrix representation
function matrix_representation(operator :: ScalarProductOperator{B, O}) :: SparseMatrixCSC{Complex{Float64}} where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    return matrix_representation(operator.op) .* operator.factor
end

# possibly recalculate the matrix representation
function recalculate!(operator :: ScalarProductOperator{B, O}, recursive::Bool=true, basis_change::Bool=true) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    # maybe recalculate recursively
    if recursive
        recalculate!(operator.op, true, basis_change)
    end
end

# set a parameter (returns (found parameter?, changed matrix?))
function set_parameter!(operator :: ScalarProductOperator{B, O}, parameter :: Symbol, value; print_result::Bool=false, recalculate::Bool=true, kwargs...) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    # check if it can be set in 1
    found_param, changed_matrix = set_parameter!(operator.op, parameter, value, print_result=print_result, recalculate=recalculate; kwargs...)
    return (found_param, changed_matrix)
end

# get a parameter (returns (found parameter?, parameter value or nothing))
function get_parameter(operator :: ScalarProductOperator{B, O}, parameter :: Symbol; print_result::Bool=false, kwargs...) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    # check if parameter can be obtained
    return get_parameter(operator.op, parameter, print_result=print_result; kwargs...)
end

# get a parameter (returns (found parameter?, parameter value or nothing))
function get_parameters(operator :: ScalarProductOperator{B, O}; kwargs...) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    # check if parameter can be obtained
    return get_parameters(operator.op; kwargs...)
end




################################################################################
#   (*) functions
################################################################################

# import of base function
import Base.(*)

# Define the (*) function for any operators
function *(f :: Number, op :: O)  where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    return ScalarProductOperator(f, op)
end
function *(op :: O, f :: Number)  where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    return ScalarProductOperator(f, op)
end
