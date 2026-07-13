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
    # for each SP state a in basis_to, precompute which SP states b in basis_from have nonzero overlap -- used to prune the MP pair loop.
    # so cross-site MP pairs are guaranteed zero and can be skipped entirely.
    sp_partners = [findall(b -> abs(sp_overlaps[a,b]) > 1e-14, 1:length(sp_from)) for a in 1:length(sp_to)]
    # reusable N x N submatrix buffer for the determinant calculation
    S = Matrix{Complex{Float64}}(undef, N, N)
    # pre-allocate reusable buffers for the inner loop
    allowed      = Int[]           # SP states in basis_from with nonzero overlap
    candidate_js = Int[]           # candidate j states before subset filter
    # triplet construction for indices and overlaps (I,J,ovl)
    Is = Int[]
    Js = Int[]
    ovls = Complex{Float64}[]
    # fill the matrix
    for i in 1:length(basis_to)
        occ_i = basis_to[i].occupation # length-N vector of SP indices
        # fill! in_allowed with false and empty! allowed
        empty!(allowed)
        # collect all SP states in basis_from that any orbital of state i could possibly overlap with
        for a in occ_i
            for b in sp_partners[a]
                if !(b in allowed)
                    push!(allowed, b)
                end
            end
        end
        # if no nonzero overlaps, skip to next iteration
        if isempty(allowed)
            continue
        end
        # collect candidate j states from lookup_sp_states for each allowed SP orbital
        empty!(candidate_js)
        for b in allowed
            for j in basis_from.lookup_sp_states[b]
                push!(candidate_js, j)
            end
        end
        unique!(candidate_js)  # deduplicate in place
        # further filter: j's entire occupation must be a subset of allowed 
        # (union above is too broad -- it includes states with orbitals outside allowed, which are guaranteed zero det)
        for j in candidate_js
            occ_j = basis_from[j].occupation  # length-N vector of SP indices
            if !all(b -> b in allowed, occ_j)
                continue
            end
            # fill using precomputed SP overlap table
            @inbounds for k in 1:N
                @inbounds for l in 1:N
                    S[k,l] = sp_overlaps[occ_i[k], occ_j[l]]
                end
            end
            # MP overlap = det(S) for fermionic states (Slater determinant overlap)
            ovl=det(S)
            if abs(ovl)> 1e-14
                push!(Is, i)
                push!(Js, j)
                push!(ovls, ovl)
            end
        end
    end
    # return the matrix
    return sparse(Is, Js, ovls, length(basis_to), length(basis_from))
end
export projector_matrix
