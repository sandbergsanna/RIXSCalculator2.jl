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

