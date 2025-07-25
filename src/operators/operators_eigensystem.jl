function eigensystem(operator :: AbstractOperator{B}; subtract_GS::Bool = false, make_hermitian::Bool=true) where {BS <: AbstractBasisState, B <: AbstractBasis{BS}}
    # build the eigensystem of the matrix
    m = matrix_representation(operator)
    if make_hermitian
        m = Hermitian((m+m')/2)
    end
    es = eigen(Matrix(m))
    evalues = Float64[real(v) for v in es.values]
    # The kth eigenvector can be obtained from the slice es.vectors[:, k]
    evectors = [es.vectors[:,k] for k in 1:length(evalues)]
    eigenlist = [(evalues[i], evectors[i]) for i in 1:length(evalues)]
    sort!(eigenlist, by=e->e[1])
    evalues = [e[1] for e in eigenlist]
    evectors = [e[2] for e in eigenlist]
    if subtract_GS
        evalues .-= minimum(evalues)
    end
    # construct dictonary
    eigsysop = Dict(
        # put in operator
        :operator => operator,
        # put in the eigenvalues
        :values => evalues,
        # put in the eigenvectors
        :vectors => evectors
    )
    # return the dictonary
    return eigsysop
end
function eigensystem(operator :: AbstractOperator{B}, basis_new::B2; subtract_GS::Bool = false, make_hermitian::Bool=true) where {BS <: AbstractBasisState, B <: AbstractBasis{BS}, BS2 <: AbstractBasisState, B2 <: AbstractBasis{BS2}}
    # build the eigensystem of the matrix
    operator = ProjectorOperator(operator, basis_new)
    m = matrix_representation(operator)
    if make_hermitian
        m = Hermitian((m+m')/2)
    end
    es = eigen(Matrix(m))
    evalues = Float64[real(v) for v in es.values]
    # The kth eigenvector can be obtained from the slice es.vectors[:, k]
    evectors = [es.vectors[:,k] for k in 1:length(evalues)]
    eigenlist = [(evalues[i], evectors[i]) for i in 1:length(evalues)]
    sort!(eigenlist, by=e->e[1])
    evalues = [e[1] for e in eigenlist]
    evectors = [e[2] for e in eigenlist]
    if subtract_GS
        evalues .-= minimum(evalues)
    end
    # construct dictonary
    eigsysop = Dict(
        # put in operator
        :operator => operator,
        # put in the eigenvalues
        :values => evalues,
        # put in the eigenvectors
        :vectors => evectors
    )
    # return the dictonary
    return eigsysop
end
export eigensystem

function toCompositeBasis(
            eigensystem :: Dict
        )
    # check if SPBasis or MPBasis
    if typeof(basis(eigensystem[:operator])) <: SPBasis
        return SPBasis([CompositeBasisState(s,basis(eigensystem[:operator])) for s in eigensystem[:vectors]])
    end
end

function printMPState(
        eigensystem :: Dict,
        index :: Integer
        ;
        subtract_GS :: Bool=false,
        cutoff :: Real = 0.01,
        digits :: Integer = 3,
        max_contributions :: Integer = 10
    )
    # get the state and the eigenenergy
    es_state  = eigensystem[:vectors][index]
    es_energy = eigensystem[:values][index]
    if subtract_GS
        es_energy -= eigensystem[:values][1]
    end
    # get the basis of the operator
    b = basis(eigensystem[:operator])
    # print the state
    printMPState(es_state, b, cutoff=cutoff, digits=digits, max_contributions=max_contributions, energy=es_energy, es_index=index)
end
export printMPState


function energies(op :: AbstractOperator)
    return sort(real.(eigvals(Matrix(matrix_representation(op)))))
end
export energies

function print_energies(op :: AbstractOperator; digits::Int64=6, subtract_GS::Bool = false)
    # get the energies
    evals = round.(energies(op), digits=digits)
    if subtract_GS
        evals .-= minimum(evals)
    end
    # save how often the energies occur
    evals_unique = unique(evals)
    evals_occurs = zeros(Int64, length(evals_unique))
    for e in evals
        for i in 1:length(evals_unique)
            if e == evals_unique[i]
                evals_occurs[i] += 1
                break
            end
        end
    end
    # print statistic
    println("$(length(evals_unique)) different energy values found:")
    for i in 1:length(evals_unique)
        println(" --> ", round(evals_unique[i], digits=6), " (x",evals_occurs[i],")")
    end
end
export print_energies

function is_on_site(
        state :: Vector{<:Number},
        site  :: Integer,
        b     :: MPBasis
        ;
        precision :: Real = 1e-3
    ) :: Bool

    # check which single particle states contain site
    state_indices = Int64[
        i for i in 1:length(b.single_particle_basis) if b.single_particle_basis[i].site == site
    ]
    # check which multi particle states contain site
    state_indices_mp = sort(unique(Int64[
        i for j in state_indices for i in b.lookup_sp_states[j]
    ]))
    # check if the state contains elements of list
    for i in state_indices_mp
        if abs(state[i]) > precision
            return true
        end
    end
    return false
end
export is_on_site

#Function that returns the indices of the eigenstates (as sorted in eigensystem) of multiplets.
#Input: eigensystem, digits_energy (optional) Output: vector of unique energies, corresponding multiplet indices
function multiplets(eigensystem :: Dict; digits_energy::Int64=6)
    # Energies and eigenvalues (already sorted)
    evals=round.(eigensystem[:values],digits=digits_energy)
    # Check unique energies
    evals_unique=unique(evals)
    # Check degeneracies and group into multiplets
    multiplet_indices=[Vector{Int64}() for _ in 1:length(evals_unique)]
    for i in 1:length(evals_unique)
        for j in 1:length(evals)
            if evals[j]==evals_unique[i]
                push!(multiplet_indices[i],j)
            end
        end
    end
    # Return unique energies and indices of corresponding multiplets
    return evals_unique,multiplet_indices
end
export multiplets

