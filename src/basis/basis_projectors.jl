# General projector as MATRICES

# projection matrix P: |basis_to> = P |basis_from>
# |basis_to><basis_from|
"""
    projector_matrix(
        basis_to   :: B1,
        basis_from :: B2
    ) :: SparseMatrixCSC{Complex{Float64}} where {BS1<:AbstractSPBasisState, BS2<:AbstractSPBasisState, B1 <: AbstractBasis{BS1}, B2 <: AbstractBasis{BS2}}

This function computes the matrix projector P, defined as:

`` \\left| basis_{to} \\right> = P \\left| basis_{from} \\right> \\longrightarrow P=\\left|basis_{to} \\right> \\left< basis_{from} \\right| ``
"""
function projector_matrix(
            basis_to   :: B1,
            basis_from :: B2
        ) :: SparseMatrixCSC{Complex{Float64}} where {BS1<:AbstractSPBasisState, BS2<:AbstractSPBasisState, B1 <: AbstractBasis{BS1}, B2 <: AbstractBasis{BS2}}
    # create new matrix
    matrix = spzeros(Complex{Float64}, length(basis_to),length(basis_from))
    # fill the matrix
    for i in 1:length(basis_to)
    for j in 1:length(basis_from)
        # fill in the element i,j (NOTE THE CONVENTION!!!)
        matrix[i,j] = overlap(basis_to[i], basis_from[j])
    end
    end
    # return the matrix
    return matrix
end
export projector_matrix


# projection matrix P: |basis_to> = P |basis_from>
# |basis_to><basis_from|
function projector_matrix(
            basis_to   :: B1,
            basis_from :: B2
        ) :: SparseMatrixCSC{Complex{Float64}} where {N,SPBS1<:AbstractSPBasisState,SPBS2<:AbstractSPBasisState, B1 <: MPBasis{N,SPBS1}, B2 <: MPBasis{N,SPBS2}}
    # get single-particle bases
    sp_to   = basis_to.single_particle_basis
    sp_from = basis_from.single_particle_basis
    # precompute the single-particle overlap matrix once
    sp_overlaps = Matrix{Complex{Float64}}(undef, length(sp_to),length(sp_from))
    for a in 1:length(sp_to)
        for b in 1:length(sp_from)
            sp_overlaps[a,b] = overlap(sp_to[a], sp_from[b])
        end
    end
    # reusable N x N submatrix buffer for the determinant calculation
    S = Matrix{Complex{Float64}}(undef, N, N)
    # create new MP matrix
    matrix = spzeros(Complex{Float64}, length(basis_to),length(basis_from))
    # fill the matrix
    for i in 1:length(basis_to)
    for j in 1:length(basis_from)
        occ_i = basis_to[i].occupation    # length-N vector of SP indices
        occ_j = basis_from[j].occupation  # length-N vector of SP indices
        # fill using precomputed SP overlap table
        @inbounds for k in 1:N
            @inbounds for l in 1:N
                S[k,l] = sp_overlaps[occ_i[k], occ_j[l]]
            end
        end
        # MP overlap = det(S) for fermionic states (Slater determinant overlap)
        v=det(S)
        if abs(v)> 1e-14
            matrix[i,j] = v
        end
    end
    end
    # return the matrix
    return matrix
end
export projector_matrix
