################################################################################
#
#   Spectrum generating functions
#
################################################################################


# Calculate the dipole amplitude when transitioning from -> D -> to
# which corresponds to <to|D|from>
function get_amplitude(
        dipole_matrix :: Matrix{Complex{Float64}},
        state_from :: Vector{<:Number},
        state_to   :: Vector{<:Number}
    ) ::Complex{Float64} where {B <: AbstractBasis}

    # build amplitude like <to|D|from>
    return dot(state_to,  dipole_matrix * state_from)
end

# Calculate the dipole amplitude when transitioning from -> D -> to
# which corresponds to <to|D|from>
function get_amplitude(
        dipole_matrix :: SparseMatrixCSC{Complex{Float64}},
        state_from :: Vector{<:Number},
        state_to   :: Vector{<:Number}
    ) ::Complex{Float64} where {B <: AbstractBasis}

    # build amplitude like <to|D|from>
    return dot(state_to,  dipole_matrix * state_from)
end

# Calculate the dipole amplitude when transitioning i_from -> D -> i_to
# which corresponds to <i_to|D|i_from>
function get_amplitude(
        eigensys        :: Dict{Symbol, Any},
        dipole_matrix :: Matrix{Complex{Float64}},
        i_from :: Int64,
        i_to   :: Int64
    ) ::Complex{Float64} where {B <: AbstractBasis}

    # build amplitude like <to|D|from>
    return get_amplitude(dipole_matrix , eigensys[:vectors][i_from], eigensys[:vectors][i_to])
end

# Calculate the dipole amplitude when transitioning i_from -> D -> i_to
# which corresponds to <i_to|D|i_from>
function get_amplitude(
        eigensys        :: Dict{Symbol, Any},
        dipole_matrix :: SparseMatrixCSC{Complex{Float64}},
        i_from :: Int64,
        i_to   :: Int64
    ) ::Complex{Float64} where {B <: AbstractBasis}

    # build amplitude like <to|D|from>
    return get_amplitude(dipole_matrix , eigensys[:vectors][i_from], eigensys[:vectors][i_to])
end
export get_amplitude


# Calculate a spectrum based on a RIXSSample object
# polarizations are in global coordinates (NOT in sample coordinates)
# basically a wrapper around getSingleSiteTransitionIntensities(sample,...)
function get_spectrum(
        eigensys        :: Dict{Symbol, Any},
        dipole_operator :: AbstractOperator{B}
        ;
        linewidth :: Real = 25.0
    ) :: Spectrum where {B <: AbstractBasis}

    # allocate the transitions
    transitions = Vector{Transition}(undef, length(eigensys[:values]))

    # get dipole matrix
    dipole_matrix=matrix_representation(dipole_operator)

    # fill all transitionss
    for i in 1:length(transitions)
        transitions[i] = Transition(
            abs( get_amplitude(eigensys, dipole_matrix, 1, i) )^2,
            linewidth,
            eigensys[:values][i] - eigensys[:values][1]
        )
    end

    # construct spectrum
    return Spectrum(transitions)
end

# get spectrum from hamiltonian with ground states in array states
function get_spectrum(
        eigensys        :: Dict{Symbol, Any},
        dipole_operator :: AbstractOperator{B},
        states          :: Vector{<:Integer}
        ;
        linewidth :: Real = 25.0
    ) :: Spectrum where {B <: AbstractBasis}
    
    # get dipole matrix
    dipole_matrix=matrix_representation(dipole_operator)

    return sum([
        begin
            # allocate the transitions
            transitions = Vector{Transition}(undef, length(eigensys[:values]))

            # fill all transitionss
            for i in 1:length(transitions)
                transitions[i] = Transition(
                    (1.0 ./ length(states)) * abs( dot(eigensys[:vectors][i],  dipole_matrix * eigensys[:vectors][j]) )^2,
                    linewidth,
                    eigensys[:values][i] - eigensys[:values][j]
                )
            end

            # construct spectrum
            Spectrum(transitions)
        end

        for j in states
    ])
end

function get_spectrum(
        eigensys        :: Dict{Symbol, Any},
        dipole_operator :: AbstractOperator{B},
        GS_energy       :: Real,
        args...
        ;
        precision :: Real = 1e-8,
        kwargs...
    ) :: Spectrum where {B <: AbstractBasis}

    # find out the states of that particular energy
    states = Int64[i for i in 1:length(eigensys[:values]) if isapprox(eigensys[:values][i], GS_energy, atol=precision)]

    # construct spectrum
    return get_spectrum(
            eigensys,
            dipole_operator,
            states,
            args...
            ;
            kwargs...
        )
end


# get spectrum from hamiltonian instead of eigensystem
function get_spectrum(
        hamiltonian     :: AbstractOperator{B},
        args...
        ;
        kwargs...
    ) :: Spectrum where {B <: AbstractBasis}

    # build eigensystem and pass on
    return get_spectrum(
        eigensystem(hamiltonian),
        args...
        ;
        kwargs...
    )
end

# export the function
export get_spectrum