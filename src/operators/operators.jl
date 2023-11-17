# using packages
using LinearAlgebra
using Combinatorics

# abstract type
include("operator_abstract_type.jl")



# single particle operators
include("operators_sp/operators_sp.jl")

# specific operators
include("specific_operators/fundamental_operators.jl")

include("specific_operators/spin_orbit_coupling.jl")

# eigensystem of an operator
include("operators_eigensystem.jl")




