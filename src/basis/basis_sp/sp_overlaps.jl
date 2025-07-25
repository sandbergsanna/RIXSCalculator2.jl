# FALLBACK DIAGONAL OVERLAP (if the diagonal overlap is not defined)
function overlap(state_1 :: SPSSBS, state_2 :: SPSSBS) :: Complex{Float64}   where {SPSSBS <: AbstractSPSSBasisState}
    # default diagonal overlap, check if the two states are identical in every field
    if state_1 == state_2
        return 1.0
    else
        return 0.0
    end
end


################################################################################
#   OVERLAP BETWEEN DIFFERENT BASIS STATES
#   NOTE THAT THE CONVENTION IS <alpha|beta> = overlap(alpha, beta)
################################################################################


# <LS|XYZ>
# = (<XYZ|LS>)^dagger
function overlap(state_1::BasisStateLS, state_2::BasisStateXYZ) :: Complex{Float64}
    # just take the conjugate of the opposite (LS -- XYZ)
    return conj(overlap(state_2, state_1))
end

# <XYZ|LS>
# with definition of t2g LS from t2g XYZ basis
# --> |+1> = -1/sqrt(2) ( |x> + i|y> )
# --> |-1> = +1/sqrt(2) ( |x> - i|y> )
# --> |0>  = |z>
# consistent with L. Ament and G. Khaliullin (Phys. Rev. B 84, 020403(R) (2011))
function overlap(state_1::BasisStateXYZ, state_2::BasisStateLS) :: Complex{Float64}
    # assert that it is an l=1,s=1/2 state (t2gs)
    @assert state_2.l == 1 && state_2.s == 1//2 "LS State is not a t2g state: $(state_2)"
    # split up depending on orbital of state 2
    if state_2.ml == 1
        # overlap with l=1 orbital
        return -1/sqrt(2.0) * (
              overlap(state_1, BasisStateXYZ(:x,state_2.ms))
            + overlap(state_1, BasisStateXYZ(:y,state_2.ms))*im
               )
    elseif state_2.ml == -1
        # overlap with l=1 orbital
        return +1/sqrt(2.0) * (
              overlap(state_1, BasisStateXYZ(:x,state_2.ms))
            - overlap(state_1, BasisStateXYZ(:y,state_2.ms))*im
               )
    elseif state_2.ml == 0
        # overlap with z orbital
        return overlap(state_1, BasisStateXYZ(:z,state_2.ms))
    else
        @error "unknown LS case, LS State appears not to be a t2g state: $(state_2)" stacktrace()
        return NaN+ NaN*im
    end
end

# <LS|J>
#
# with definition of J basis
# --> |1/2,-1/2> = sqrt(1/3) | 0,-> - sqrt(2/3) |-1,+>
# --> |1/2,+1/2> = sqrt(1/3) | 0,+> - sqrt(2/3) |+1,->
# --> |3/2,-3/2> = |-1,->
# --> |3/2,-1/2> = sqrt(2/3) |0,-> + sqrt(1/3) |-1,+>
# --> |3/2,+1/2> = sqrt(2/3) |0,+> + sqrt(1/3) |+1,->
# --> |3/2,+3/2> = |+1,+>

function overlap(state_1::Union{BasisStateLS, BasisStateXYZ}, state_2::BasisStateJ) :: Complex{Float64}
    # split up depending on numbers of state 2
    if state_2.j == 1//2
        # j=1//2, depending on mjp
        if state_2.mj == -1//2
            # overlap with \1//2, -1//2> state
            return (
                  sqrt(1.0/3.0) * overlap(state_1, BasisStateLS(1, 0, 1//2,-1//2))
                - sqrt(2.0/3.0) * overlap(state_1, BasisStateLS(1,-1, 1//2,+1//2))
                   )
        elseif state_2.mj == +1//2
            # overlap with \1//2, +1//2> state
            return (
                 sqrt(1.0/3.0) * overlap(state_1, BasisStateLS(1, 0, 1//2, 1//2))
               - sqrt(2.0/3.0) * overlap(state_1, BasisStateLS(1, 1, 1//2,-1//2))
                  )
        else
            @error "wrong value of mj = $(state_2.mj)" stacktrace()
            return NaN+ NaN*im
        end
    elseif state_2.j == 3//2
        # j=3//2, depending on mjp
        if state_2.mj == -3//2
            # overlap with \3//2, -3//2> state
            return overlap(state_1, BasisStateLS(1,-1, 1//2,-1//2))
        elseif state_2.mj == -1//2
            # overlap with \3//2, -1//2> state
            return (
                  sqrt(1.0/3.0) * overlap(state_1, BasisStateLS(1,-1, 1//2,+1//2))
                + sqrt(2.0/3.0) * overlap(state_1, BasisStateLS(1, 0, 1//2,-1//2))
                  )
        elseif state_2.mj == 1//2
            # overlap with \3//2, 1//2> state
            return (
                  sqrt(1.0/3.0) * overlap(state_1, BasisStateLS(1, 1, 1//2,-1//2))
                + sqrt(2.0/3.0) * overlap(state_1, BasisStateLS(1, 0, 1//2,+1//2))
                )
        elseif state_2.mj == 3//2
            # overlap with \3//2, +3//2> state
            return overlap(state_1, BasisStateLS(1, 1, 1//2, 1//2))
        else
            @error "wrong value of mj = $(state_2.mj)" stacktrace()
            return NaN+ NaN*im
        end
    else
        @error "uknown value of j = $(state_2.j)" stacktrace()
        return NaN+ NaN*im
    end
end

# <J|LS>
# = (<LS|J>)^dagger
function overlap(state_1::BasisStateJ, state_2::Union{BasisStateLS, BasisStateXYZ}) :: Complex{Float64}
    # just take the conjugate of the opposite (LS -- J)
    return conj(overlap(state_2, state_1))
end




# <alpha|A1G>
# with definition of electron A1G basis
#   --> |a1g> = 1/sqrt(3) ( |x> + |y> + |z> )
#   --> |eg+> = 1/sqrt(3) ( e^(-im 2pi/3)|x> + e^(+im 2pi/3)|y> + |z> )
#   --> |eg-> =-1/sqrt(3) ( e^(+im 2pi/3)|x> + e^(-im 2pi/3)|y> + |z> )
# see Khomskii Book, Equation (3.13) and transform to hole picture (cc)
#   --> |eg1> = -im/sqrt(2) ( |eg+> + |eg-> )
#   --> |eg2> =   1/sqrt(2) ( |eg+> - |eg-> )
# see Khomskii Book, Equation (3.14, 3.15)
function overlap(state_1::Union{BasisStateXYZ,BasisStateJ,BasisStateLS,BasisStateA1G}, state_2::BasisStateA1G) :: Complex{Float64}
    # split up depending on numbers of state 1
    if state_2.orbital == :a1g
        # overlap with x orbital
        return 1/sqrt(3.0) * (
              overlap(state_1, BasisStateXYZ(:x, state_2.ms))
            + overlap(state_1, BasisStateXYZ(:y, state_2.ms))
            + overlap(state_1, BasisStateXYZ(:z, state_2.ms))
               )
    elseif state_2.orbital == :egp
        # overlap with x orbital
        return 1/sqrt(3.0) * (
              overlap(state_1, BasisStateXYZ(:x, state_2.ms)) * exp( im * 2pi * (1.0/3.0))
            + overlap(state_1, BasisStateXYZ(:y, state_2.ms)) * exp( im * 2pi * (2.0/3.0))
            + overlap(state_1, BasisStateXYZ(:z, state_2.ms)) * exp( im * 2pi * (3.0/3.0))
               )
    elseif state_2.orbital == :egm
        # overlap with x orbital
        return -1/sqrt(3.0) * (
              overlap(state_1, BasisStateXYZ(:x, state_2.ms)) * exp(-im * 2pi * (1.0/3.0))
            + overlap(state_1, BasisStateXYZ(:y, state_2.ms)) * exp(-im * 2pi * (2.0/3.0))
            + overlap(state_1, BasisStateXYZ(:z, state_2.ms)) * exp(-im * 2pi * (3.0/3.0))
               )
    elseif state_2.orbital == :eg1
        # overlap with x orbital
        return -im/sqrt(2.0) * (
              overlap(state_1, BasisStateA1G(:egp, state_2.ms))
            + overlap(state_1, BasisStateA1G(:egm, state_2.ms))
               )
    elseif state_2.orbital == :eg2
        # overlap with x orbital
        return 1/sqrt(2.0) * (
              overlap(state_1, BasisStateA1G(:egp, state_2.ms))
            - overlap(state_1, BasisStateA1G(:egm, state_2.ms))
               )
    else
        @error "unknown a1g orbital $(state_2.orbital)" stacktrace()
        return NaN+ NaN*im
    end
end

# <A1G|alpha>
# = (<alpha|A1G>)^dagger
function overlap(state_1::BasisStateA1G, state_2::Union{BasisStateXYZ,BasisStateJ,BasisStateLS}) :: Complex{Float64}
    # just take the conjugate of the opposite (XYZ -- A1G)
    return conj(overlap(state_2, state_1))
end



###################################################################################
# Definition of overlap for SINGLE PARTICLE SINGLE SITE COMPOSITE Basis states #
####################################################################################

# Definition of overlap for SINGLE PARTICLE SINGLE SITE and SINGLE PARTICLE SINGLE SITE  COMPOSITE Basis states
# <alpha|COMP> = prefactor_j <alpha|basisstate_j>
function overlap(state_1 :: SPSSP, state_2 :: SPSSCompositeBasisState{<:SPBasis{SPSS}}) :: Complex{Float64} where {SPSS<:AbstractSPSSBasisState, SPSSP<:AbstractSPSSBasisState}
    # write down the overlap
    ovl = 0.0+0.0im
    # overlap with all basis states
    for bsi in eachindex(state_2.basis)
        ovl +=  state_2.prefactors[bsi] * overlap(state_1, state_2.basis[bsi])
    end
    # return the overlap
    return ovl
end

# Definition of overlap for SINGLE PARTICLE SINGLE SITE and SINGLE PARTICLE SINGLE SITE  COMPOSITE Basis states
# <COMP|alpha> = (<alpha|COMP>)^dagger
function overlap(state_1 :: SPSSCompositeBasisState{<:SPBasis{SPSS}}, state_2 :: SPSSP) :: Complex{Float64} where {SPSS<:AbstractSPSSBasisState, SPSSP<:AbstractSPSSBasisState}
    # return the conjugate
    return conj(overlap(state_2, state_1))
end



# Definition of overlap for SINGLE PARTICLE SINGLE SITE COMPOSITE Basis states
# <COMP1|COMP2> = sum_j prefactor2_j <COMP1|basisstate2_j>
function overlap(state_1 :: SPSSCompositeBasisState{<:SPBasis{SPSS}}, state_2 :: SPSSCompositeBasisState{<:SPBasis{SPSSP}}) :: Complex{Float64} where {SPSS<:AbstractSPSSBasisState, SPSSP<:AbstractSPSSBasisState}
    # write down the overlap
    ovl = 0.0+0.0im
    # overlap with all basis states
    for bsi in eachindex(state_2.basis)
        ovl +=  state_2.prefactors[bsi] * overlap(state_1, state_2.basis[bsi])
    end
    # return the overlap
    return ovl
end

#####################################################################
# Definition of overlap for SINGLE PARTICLE MULTI SITE Basis states #
#####################################################################

# Definition of overlap for SINGLE PARTICLE MULTI SITE Basis states
# <i,alpha|j,beta> = delta(i,j) * <alpha|beta>
function overlap(state_1 :: SPMSBasisState{BS1}, state_2 :: SPMSBasisState{BS2}) :: Complex{Float64} where {BS1<:AbstractSPSSBasisState, BS2<:AbstractSPSSBasisState}
    # check if the states are on the same site,
    if state_1.site == state_2.site
        # if they are, return their overlap
        return overlap(state_1.state, state_2.state)
    else
        # if not return 0
        return 0.0
    end
end

# Definition of overlap for SINGLE PARTICLE MULTI SITE Basis states
# <alpha|COMP> = prefactor_j <alpha|basisstate_j>
function overlap(state_1 :: SPMSP, state_2 :: SPMSCompositeBasisState{<:SPBasis{<:SPMSBasisState{SPSS}}}) :: Complex{Float64} where {SPSS<:AbstractSPSSBasisState, SPSSP<:AbstractSPSSBasisState, SPMSP<:SPMSBasisState{SPSSP}}
    # write down the overlap
    ovl = 0.0+0.0im
    # overlap with all basis states
    for bsi in eachindex(state_2.basis)
        ovl +=  state_2.prefactors[bsi] * overlap(state_1, state_2.basis[bsi])
    end
    # return the overlap
    return ovl
end

# Definition of overlap for SINGLE PARTICLE MULTI SITE Basis states
# <COMP|alpha> = (<alpha|COMP>)^dagger
function overlap(state_1 :: SPMSCompositeBasisState{<:SPBasis{<:SPMSBasisState{SPSS}}}, state_2 :: SPMSP) :: Complex{Float64} where {SPSS<:AbstractSPSSBasisState, SPSSP<:AbstractSPSSBasisState, SPMSP<:SPMSBasisState{SPSSP}}
    # return the conjugate
    return conj(overlap(state_2, state_1))
end

# Definition of overlap for SINGLE PARTICLE MULTI SITE Basis states
# <COMP1|COMP2> = sum_j prefactor2_j <COMP1|basisstate2_j>
function overlap(state_1 :: SPMSCompositeBasisState{<:SPBasis{<:SPMSBasisState{SPSS}}}, state_2 :: SPMSCompositeBasisState{<:SPBasis{<:SPMSBasisState{SPSSP}}}) :: Complex{Float64} where {SPSS<:AbstractSPSSBasisState, SPSSP<:AbstractSPSSBasisState}
    # write down the overlap
    ovl = 0.0+0.0im
    # overlap with all basis states
    for bsi in eachindex(state_2.basis)
        ovl +=  state_2.prefactors[bsi] * overlap(state_1, state_2.basis[bsi])
    end
    # return the overlap
    return ovl
end
