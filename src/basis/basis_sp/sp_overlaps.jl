# FALLBACK DIAGONAL OVERLAP (if the diagonal overlap is not defined)
function overlap(state_1 :: SPSSBS, state_2 :: SPSSBS) :: Complex{Float64}   where {SPSSBS <: AbstractSPSSBasisState}
    # default diagonal overlap, check if the two states are identical in every field
    if state_1 == state_2
        return 1.0
    else
        return 0.0
    end
end




#####################################################################
# Definition of overlap for SINGLE PARTICLE MULTI SITE Basis states #
#####################################################################


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

# Definition of overlap for SINGLE PARTICLE MULTI SITE Basis states
# <COMP|alpha> = (<alpha|COMP>)^dagger
function overlap(state_1 :: SPSSCompositeBasisState{<:SPBasis{SPSS}}, state_2 :: SPSSP) :: Complex{Float64} where {SPSS<:AbstractSPSSBasisState, SPSSP<:AbstractSPSSBasisState}
    # return the conjugate
    return conj(overlap(state_2, state_1))
end



# Definition of overlap for SINGLE PARTICLE MULTI SITE Basis states
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
