# obtaining a spectrum
function get_spectrum(
            ls :: LabSystem,
            args...
            ;
            kwargs...
        ) :: Spectrum where {B <: AbstractBasis}

    # construct spectrum
    return get_spectrum(ls.eigensys, ls.dipole_hor, args...; kwargs...) + get_spectrum(ls.eigensys, ls.dipole_ver, args...; kwargs...)
end

"""
   dq_dependence_multiplet(
        lab::LabSystem,
        dq_values::Vector{<:Real},
        q_beam::Real,
        to_multiplet::Int64
    ) 
Function that calculates intensities vs transferred momentum dq (from gs multiplet to a given multiplet). 
q_beam is the magnitude of q_in and q_out.
"""
function dq_dependence_multiplet(lab::LabSystem,dq_values::Vector{<:Real},q_beam::Real,to_multiplet::Int64)
    # get multiplets
    energy_values,multiplet_indices=multiplets(lab.eigensys)
    # initalize intensities
    intensities=zeros(length(dq_values))
    # Iterate over all values for dQ
    for i in 1:length(dq_values)
        set_dQ!(lab,dq_values[i],q_beam)
        recalculate_dipole_operators!(lab)
        # get dipole matrices
        dipole_matrix_hor=matrix_representation(lab.dipole_hor)
        dipole_matrix_ver=matrix_representation(lab.dipole_ver)
        # Iterate over states in excited multiplet
        for j in multiplet_indices[to_multiplet]
            # Iterate over states in gs multiplet
            for k in multiplet_indices[1]
                intensities[i]+=abs(get_amplitude(lab.eigensys,dipole_matrix_hor,k,j)+get_amplitude(lab.eigensys,dipole_matrix_ver,k,j))^2
            end
        end
    end
    return intensities
end
export dq_dependence_multiplet

"""
    theta_dependence_multiplet(
        lab::LabSystem,
        theta_values::Vector{<:Real}, 
        twotheta :: Real, 
        q_beam :: Real, 
        to_multiplet::Int64
    )
Function that calculates intensities vs theta (from gs multiplet to a given multiplet). 
q_beam is the magnitude of q_in and q_out.
"""
function theta_dependence_multiplet(lab::LabSystem,theta_values::Vector{<:Real}, twotheta :: Real, q_beam :: Real, to_multiplet::Int64)
    # get multiplets
    energy_values,multiplet_indices=multiplets(lab.eigensys)
    # initalize intensities
    intensities=zeros(length(theta_values))
    # Iterate over all values for theta
    for i in 1:length(theta_values)
         # set scattering angles
        set_scattering_angles_deg!(lab, theta_values[i],twotheta, q_beam)
        recalculate_dipole_operators!(lab)
        # get dipole matrices
        dipole_matrix_hor=matrix_representation(lab.dipole_hor)
        dipole_matrix_ver=matrix_representation(lab.dipole_ver)
        # Iterate over states in excited multiplet
        for j in multiplet_indices[to_multiplet]
            # Iterate over states in gs multiplet
            for k in multiplet_indices[1]
                intensities[i]+=abs(get_amplitude(lab.eigensys,dipole_matrix_hor,k,j)+get_amplitude(lab.eigensys,dipole_matrix_ver,k,j))^2
            end
        end
    end
    return intensities
end

"""
    theta_dependence_multiplet(
        lab::LabSystem,
        theta_values::Vector{<:Real}, 
        twotheta_values :: Vector{<:Real}, 
        q_beam :: Real, 
        to_multiplet::Int64
    )
Function that calculates intensities vs theta (from gs multiplet to a given multiplet). 
q_beam is the magnitude of q_in and q_out.
"""
function theta_dependence_multiplet(lab::LabSystem,theta_values::Vector{<:Real}, twotheta_values :: Vector{<:Real}, q_beam :: Real, to_multiplet::Int64)
    # get multiplets
    energy_values,multiplet_indices=multiplets(lab.eigensys)
    # initalize intensities
    intensities=zeros(length(theta_values))
    # Iterate over all values for theta and twotheta values
    for i in 1:length(theta_values)
         # set scattering angles
        set_scattering_angles_deg!(lab, theta_values[i],twotheta_values[i], q_beam)
        recalculate_dipole_operators!(lab)
        # get dipole matrices
        dipole_matrix_hor=matrix_representation(lab.dipole_hor)
        dipole_matrix_ver=matrix_representation(lab.dipole_ver)
        # Iterate over states in excited multiplet
        for j in multiplet_indices[to_multiplet]
            # Iterate over states in gs multiplet
            for k in multiplet_indices[1]
                intensities[i]+=abs(get_amplitude(lab.eigensys,dipole_matrix_hor,k,j)+get_amplitude(lab.eigensys,dipole_matrix_ver,k,j))^2
            end
        end
    end
    return intensities
end

export theta_dependence_multiplet