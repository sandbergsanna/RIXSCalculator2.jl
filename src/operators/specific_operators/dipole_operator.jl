# Type Definition of Dipole Operator
"""
    mutable struct DipoleOperator <: AbstractSPMSOperator{SPBasis{SPMSBasisState{BasisStateXYZ}}}

`DipoleOperator` refers to the dipole operator.

# Fields

- `basis :: SPBasis{SPMSBasisState{BasisStateXYZ}}`, a single particle basis ;

- `edge :: Int64`, the current edge;

- `spin_quantization:: CoordinateFrame`, the spin quantization frame;

- `site :: Int64`, the site;  

- `position :: Vector{Float64}`, its position;

- `eps_in :: Vector{Float64}`, the ingoing polarization;

- `q_in :: Vector{Float64}`, the ingoing momentum;

- `eps_out :: Vector{Float64}`, the outgoing polarization;

- `q_out :: Vector{Float64}`, the ingoing momentum;

"""
mutable struct DipoleOperator <: AbstractSPMSOperator{SPBasis{SPMSBasisState{BasisStateXYZ}}}

    # the inner basis (XYZ type basis)
    basis :: SPBasis{SPMSBasisState{BasisStateXYZ}}

    # the current edge
    edge :: Int64

    # spin quantization frame
    spin_quantization :: CoordinateFrame

    # the site index
    site :: Int64

    # the site position
    position :: Vector{Float64}

    # the ingoing polarization
    eps_in  :: Vector{Float64}
    # the ingoing momentum
    q_in    :: Vector{Float64}

    # the outgoing polarization
    eps_out :: Vector{Float64}
    # the outgoing momentum
    q_out   :: Vector{Float64}

end

"""
    DipoleOperator(
                ;
                edge::Int64=3,
                site::Int64=1,
                site_position::Vector{<:Real}=[0,0,0],
                eps_in  ::Vector{<:Real}=[1,0,0],
                eps_out ::Vector{<:Real}=[1,0,0],
                q_in    ::Vector{<:Real}=[0,0,0],
                q_out   ::Vector{<:Real}=[0,0,0]
            ) :: DipoleOperator

Creates an object of type `DipoleOperator` with the input values.
"""
function DipoleOperator(
            ;
            edge::Int64=3,
            site::Int64=1,
            site_position::Vector{<:Real}=[0,0,0],
            eps_in  ::Vector{<:Real}=[1,0,0],
            eps_out ::Vector{<:Real}=[1,0,0],
            q_in    ::Vector{<:Real}=[0,0,0],
            q_out   ::Vector{<:Real}=[0,0,0]
        ) :: DipoleOperator

    # create a new object
    op = DipoleOperator(
        SPBasis{SPMSBasisState{BasisStateXYZ}}([ SPMSBasisState{BasisStateXYZ}(state, site) for state in states(getT2GBasisXYZ()) ]),
        edge,
        CoordinateFrame(),
        site,
        site_position,
        eps_in,
        q_in,
        eps_out,
        q_out
    )

    # return the object
    return op
end

function DipoleOperator(
        basis :: SPBasis{SPMSBasisState{SPSS}}
        ;
        kwargs...
    ) where {SPSS <: AbstractSPSSBasisState}

    # return the projection to that basis
    return ProjectorOperator(DipoleOperator(;kwargs...), basis)
end

function DipoleOperator(
        basis :: MPB
        ;
        kwargs...
    ) where {N, SPSS <: AbstractSPSSBasisState, SPBS <: SPMSBasisState{SPSS}, MPB <: MPBasis{N,SPBS}}

    # return the 1p generalization of projection to that basis
    return MPGeneralizedSPOperator(basis, ProjectorOperator(DipoleOperator(;kwargs...), basis.single_particle_basis))
end
export DipoleOperator







################################################################################
#   Interface functions
################################################################################

# show the operator
import Base.show
function Base.show(io::IO, op :: DipoleOperator)
    if haskey(io, :compact)
        print(io, "L=$(op.edge) edge Dipole operator @site("*string(op.site)*") ")
    else
        print(io, "Dipole operator on site "*string(op.site)*" using L=$(op.edge) edge\n")
    end
end

# obtain the current basis
function basis(op :: DipoleOperator)
    return op.basis
end

# calculate the matrix representation
function matrix_representation(op :: DipoleOperator) :: SparseMatrixCSC{Complex{Float64}}
    # calculate the bases
    # 1) calculate the p basis
    # define the spin orbit operator in the relevant quantization axis
    SO_operator = SpinOrbitOperator(getT2GBasisLS(), 10.0)
    set_parameter!(SO_operator, :spin_quantization, op.spin_quantization)
    # generate spin orbit eigensystem
    so_es  = eigensystem(SO_operator)
    # generate the p basis out of this
    p_basis = deepcopy(toCompositeBasis(so_es))
    # set the core basis accordingly
    if op.edge == 2
        # p = 1/2 wanted
        basis_core = SPBasis{SPMSBasisState{SPSSCompositeBasisState{SPBasis{BasisStateLS}}}}(
            [ SPMSBasisState{SPSSCompositeBasisState{SPBasis{BasisStateLS}}}(state, op.site) for state in states(p_basis)[1:2]]
        )
    elseif op.edge == 3
        # p = 3/2 wanted
        basis_core = SPBasis{SPMSBasisState{SPSSCompositeBasisState{SPBasis{BasisStateLS}}}}(
            [ SPMSBasisState{SPSSCompositeBasisState{SPBasis{BasisStateLS}}}(state, op.site) for state in states(p_basis)[3:6]]
        )
    elseif op.edge == -1
        # no p selection wanted
        basis_core = SPBasis{SPMSBasisState{SPSSCompositeBasisState{SPBasis{BasisStateLS}}}}(
            [ SPMSBasisState{SPSSCompositeBasisState{SPBasis{BasisStateLS}}}(state, op.site) for state in states(p_basis) ]
        )
    else
        error("Unknown edge: ", op.edge)
    end

     # 2) recalculate the t2g basis
    op.basis = SPBasis{SPMSBasisState{BasisStateXYZ}}(
        [ SPMSBasisState{BasisStateXYZ}(state, op.site) for state in states(getT2GBasisXYZ()) ]
    )

    # calculate the internal matrices
    # 1) recalculate the core hole projection
    core_hole_projection =  projector_matrix(op.basis, basis_core) * projector_matrix(basis_core, op.basis)

    # 2) calculate the D matrices
    D_x,D_y,D_z=D_matrices()

    # normalize polarizations
    op.eps_in  = op.eps_in  ./ norm(op.eps_in)
    op.eps_out = op.eps_out ./ norm(op.eps_out)

    # calculate the matrix representation

    # compose D matrices
    D_in  = ((D_x .* op.eps_in[1])   .+
             (D_y .* op.eps_in[2])   .+
             (D_z .* op.eps_in[3]))  .* exp(-im*dot(op.position, op.q_in))
    D_out = ((D_x .* op.eps_out[1])  .+
             (D_y .* op.eps_out[2])  .+
             (D_z .* op.eps_out[3])) .* exp(-im*dot(op.position, op.q_out))

    # compose the complete string
    matrix_rep = D_out' * core_hole_projection * D_in

    # return nothing
    return matrix_rep
end

# calculate D matrices 
function D_matrices() :: Tuple{SparseMatrixCSC{Complex{Float64}},SparseMatrixCSC{Complex{Float64}},SparseMatrixCSC{Complex{Float64}}}
    # calculate the D matrices
    D_x = spzeros(Complex{Float64}, 6,6)
    D_y = spzeros(Complex{Float64}, 6,6)
    D_z = spzeros(Complex{Float64}, 6,6)
    # fill the matrices
    D_x[3,5] = 1    # <y|x|z> down
    D_x[5,3] = 1    # <z|x|y> down
    D_x[4,6] = 1    # <y|x|z> up
    D_x[6,4] = 1    # <z|x|y> up
    D_y[1,5] = 1    # <x|y|z> down
    D_y[5,1] = 1    # <z|y|x> down
    D_y[2,6] = 1    # <x|y|z> up
    D_y[6,2] = 1    # <z|y|x> up
    D_z[1,3] = 1    # <x|z|y> down
    D_z[3,1] = 1    # <y|z|x> down
    D_z[2,4] = 1    # <x|z|y> up
    D_z[4,2] = 1    # <y|z|x> up

    # return D matrices
    return D_x,D_y,D_z
end