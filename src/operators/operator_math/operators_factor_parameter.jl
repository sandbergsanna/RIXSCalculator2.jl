# Type Definition
"""
    mutable struct SettableScalarProductOperator{B, O<:AbstractOperator{B}} <: AbstractOperator{B}

This object defines the (settable) scalar product operator.

# Fields

- `basis :: B`, the basis;
- `factor :: Complex{Float64}`, `label  :: Symbol`, the scalar factor and its label;
- `op :: O`, the contained operator

"""
mutable struct SettableScalarProductOperator{B, O<:AbstractOperator{B}} <: AbstractOperator{B}
    # the basis
    basis :: B
    # the scalar factor
    factor :: Complex{Float64}
    label  :: Symbol
    # the contained operator
    op :: O

    # Custom constructor (without explicit matrix rep)
    function SettableScalarProductOperator(label::Symbol, factor::Number, op::O) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
        # construct new operator
        operator = new{B, O}(basis(op), factor, label, op)
        # return the operator
        return operator
    end
    # Custom constructor (without explicit matrix rep)
    function SettableScalarProductOperator(label::Symbol, op::O) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
        # construct new operator
        operator = new{B, O}(basis(op), 1.0, label, op)
        # return the operator
        return operator
    end
end

# export operator type
export SettableScalarProductOperator

import Base.show
function Base.show(io::IO, op::SettableScalarProductOperator{B, O}) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    if haskey(io, :compact)
        print(io, "($(op.label)=$(op.factor)) * {")
        show(io, op.op)
        print(io, "}")
    else
        print(io, "($(op.label)=$(op.factor)) times ")
        show(io, op.op)
    end
end




################################################################################
#   Interface functions
################################################################################

# obtain the current basis
function basis(operator :: SettableScalarProductOperator{B, O}) :: B where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    return operator.basis
end

# calculate the matrix representation
function matrix_representation(operator :: SettableScalarProductOperator{B, O}) :: SparseMatrixCSC{Complex{Float64}} where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    return matrix_representation(operator.op) .* operator.factor
end

# possibly recalculate the matrix representation
function recalculate!(operator :: SettableScalarProductOperator{B, O}, basis_change::Bool=true) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    # recalculate conatined opeerator
    recalculate!(operator.op, basis_change)
end

# set a parameter (returns (found parameter?, changed matrix?))
function set_parameter!(operator :: SettableScalarProductOperator{B, O}, parameter :: Symbol, value; print_result::Bool=false, recalculate::Bool=true, kwargs...) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    # check if it can be set in 1
    found_param, changed_matrix = set_parameter!(operator.op, parameter, value, print_result=print_result, recalculate=recalculate; kwargs...)
    if parameter == operator.label
        found_param_i = true
        operator.factor = value
    else
        found_param_i = false
    end
    return (found_param || found_param_i, found_param_i || changed_matrix)
end

# get a parameter (returns (found parameter?, parameter value or nothing))
function get_parameter(operator :: SettableScalarProductOperator{B, O}, parameter :: Symbol; print_result::Bool=false, kwargs...) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    # check if parameter can be obtained
    param_1 = (parameter == operator.label) ? operator.factor : nothing
    param_2 = get_parameter(operator.op, parameter, print_result=print_result; kwargs...)
    # give warning if paramter found in both
    if param_1 != nothing && param_2 != nothing
        if print_result
            println("WARNING, parameter :$(parameter) found both in sub operator and factor operator. Returning value of factor operator")
        end
        return param_1
    elseif param_1 != nothing
        if print_result
            println("Parameter :$(parameter) found in factor operator, value $(param_1)")
        end
        return param_1
    elseif param_2 != nothing
        if print_result
            println("Parameter :$(parameter) found in sub operator, value $(param_2)")
        end
        return param_2
    else
        if print_result
            println("Parameter :$(parameter) not found in factor operator or sub operator")
        end
    end
end

# get a parameter (returns (found parameter?, parameter value or nothing))
function get_parameters(operator :: SettableScalarProductOperator{B, O}; kwargs...) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    # check if parameter can be obtained
    l = get_parameters(operator.op; kwargs...)
    push!(l, operator.label)
    return unique(l)
end




################################################################################
#   (*) functions
################################################################################

# import of base function
import Base.(*)

# Define the (*) function for any operators
function *(f :: Symbol, op :: O)  where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    return SettableScalarProductOperator(f, op)
end
function *(op :: O, f :: Symbol)  where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    return SettableScalarProductOperator(f, op)
end


# Define the (*) function for any operators
function *(f :: Pair{Symbol,<:Number}, op :: O)  where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    return SettableScalarProductOperator(f[1], f[2], op)
end
function *(op :: O, f :: Pair{Symbol,<:Number})  where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    return SettableScalarProductOperator(f[1], f[2], op)
end
