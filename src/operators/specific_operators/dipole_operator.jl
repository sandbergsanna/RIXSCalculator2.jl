# Type Definition of Dipole Operator
"""
    mutable struct DipoleOperator <: AbstractSPMSOperator{SPBasis{SPMSBasisState{BasisStateXYZ}}}

`DipoleOperator` refers to the dipole operator.

# Fields

- `basis :: SPBasis{SPMSBasisState{BasisStateXYZ}}`, a single particle basis ;

- `edge :: Int64`, the current edge;

- `basis_core :: SPBasis{SPMSBasisState{SPSSCompositeBasisState{SPBasis{BasisStateLS}}}}`, the core hole basis  (p orbitals or core with j=1/2 or j=3/2);

- `spin_quantization:: CoordinateFrame`, the spin quantization frame;

- `core_hole_projection :: Matrix{Complex{Float64}}`, the core-hole projection;

- `site :: Int64`, the site;  

- `position :: Vector{Float64}`, its position;

- `eps_in :: Vector{Float64}`, the ingoing polarization;

- `q_in :: Vector{Float64}`, the ingoing momentum;

- `eps_out :: Vector{Float64}`, the outgoing polarization;

- `q_out :: Vector{Float64}`, the ingoing momentum;

- `D_x`, `D_y`, `D_z`, dipole matrices all of type `Matrix{Complex{Float64}}`;

- `matrix_rep :: Matrix{Complex{Float64}}`, the current matrix representation without prefactor.

"""
mutable struct DipoleOperator <: AbstractSPMSOperator{SPBasis{SPMSBasisState{BasisStateXYZ}}}

    # the inner basis (XYZ type basis)
    basis :: SPBasis{SPMSBasisState{BasisStateXYZ}}

    # the current edge
    edge :: Int64
    # the core hole basis (p orbitals or core with j=1/2 or j=3/2)
    basis_core :: SPBasis{SPMSBasisState{SPSSCompositeBasisState{SPBasis{BasisStateLS}}}}

    # spin quantization frame
    spin_quantization :: CoordinateFrame

    # the core hole projection (L=2 / L=3)
    core_hole_projection :: Matrix{Complex{Float64}}

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

    # dipole matrices
    D_x :: Matrix{Complex{Float64}}
    D_y :: Matrix{Complex{Float64}}
    D_z :: Matrix{Complex{Float64}}

    # the current matrix representation (without prefactor)
    matrix_rep :: Matrix{Complex{Float64}}

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
    # get D matrices
    D_x,D_y,D_z=D_matrices()
    # create a new object
    op = DipoleOperator(
        SPBasis{SPMSBasisState{BasisStateXYZ}}([ SPMSBasisState{BasisStateXYZ}(state, site) for state in states(getT2GBasisXYZ()) ]),
        edge,
        SPBasis{SPMSBasisState{SPSSCompositeBasisState{SPBasis{BasisStateLS}}}}([]),
        CoordinateFrame(),
        zeros(Complex{Float64},2,2),
        site,
        site_position,
        eps_in,
        q_in,
        eps_out,
        q_out,
        D_x,
        D_y,
        D_z,
        zeros(Complex{Float64},6,6)
    )

    # recalculate everything
    recalculate!(op,true)

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

# obtain the matrix representation
function matrix_representation(op :: DipoleOperator) :: Matrix{Complex{Float64}}
    return op.matrix_rep
end

# recalculate the matrix
function recalculate!(op :: DipoleOperator, basis_change::Bool=true)
    # recalculate the bases
    if basis_change

        # 1) recalculate the p basis
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
            op.basis_core = SPBasis{SPMSBasisState{SPSSCompositeBasisState{SPBasis{BasisStateLS}}}}(
                [ SPMSBasisState{SPSSCompositeBasisState{SPBasis{BasisStateLS}}}(state, op.site) for state in states(p_basis)[1:2]]
            )
        elseif op.edge == 3
            # p = 3/2 wanted
            op.basis_core = SPBasis{SPMSBasisState{SPSSCompositeBasisState{SPBasis{BasisStateLS}}}}(
                [ SPMSBasisState{SPSSCompositeBasisState{SPBasis{BasisStateLS}}}(state, op.site) for state in states(p_basis)[3:6]]
            )
        elseif op.edge == -1
            # no p selection wanted
            op.basis_core = SPBasis{SPMSBasisState{SPSSCompositeBasisState{SPBasis{BasisStateLS}}}}(
                [ SPMSBasisState{SPSSCompositeBasisState{SPBasis{BasisStateLS}}}(state, op.site) for state in states(p_basis) ]
            )
        else
            error("Unknown edge: ", op.edge)
        end

        # 2) recalculate the t2g basis
        op.basis = SPBasis{SPMSBasisState{BasisStateXYZ}}(
            [ SPMSBasisState{BasisStateXYZ}(state, op.site) for state in states(getT2GBasisXYZ()) ]
        )

        # 3) recalculate the core hole projection
        op.core_hole_projection =  projector_matrix(op.basis, op.basis_core) * projector_matrix(op.basis_core, op.basis)

    end

    # normalize polarizations
    op.eps_in  = op.eps_in  ./ norm(op.eps_in)
    op.eps_out = op.eps_out ./ norm(op.eps_out)

    # recalculate the matrix representation
    # compose D matrices
    D_in  = ((op.D_x .* op.eps_in[1])   .+
             (op.D_y .* op.eps_in[2])   .+
             (op.D_z .* op.eps_in[3]))  .* exp(-im*dot(op.position, op.q_in))
    D_out = ((op.D_x .* op.eps_out[1])  .+
             (op.D_y .* op.eps_out[2])  .+
             (op.D_z .* op.eps_out[3])) .* exp(-im*dot(op.position, op.q_out))

    # compose the complete string
    op.matrix_rep = D_out' * op.core_hole_projection * D_in

    # return nothing
    return nothing
end

#####################################################################
#                    D matrices
####################################################################

# calculate D matrices 
function D_matrices() :: Tuple{Matrix{Complex{Float64}},Matrix{Complex{Float64}},Matrix{Complex{Float64}}}
    # calculate the D matrices
    D_x = zeros(Complex{Float64}, 6,6)
    D_y = zeros(Complex{Float64}, 6,6)
    D_z = zeros(Complex{Float64}, 6,6)
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