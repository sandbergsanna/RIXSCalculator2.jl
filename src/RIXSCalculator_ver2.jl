################################################################################
#
#   GENERAL MODULE FOR RIXS CALCULATIONS
#   AUTHOR Jan Attig
#          Luca Peterlini
#          Anna Sandberg (2.0)
#   JULIA v.1+
#
################################################################################



# start of module
module RIXSCalculator_ver2


    # GENERAL USINGS
    using LinearAlgebra



    # INCLUDE EVERTHING FROM SUBFILES
    include("basis/basis.jl")
    include("coordinate_frames/coordinate_frames.jl")
    include("operators/operators.jl")
    


# end of module
end