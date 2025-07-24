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

#Function - Calculate intensities vs transferred momentum dq (from gs multiplet to given multiplet). 
# Input:LabSystem, dq_values, q_beam, to_multiplet
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

#Function - Calculate intensities vs theta (from gs multiplet to a given multiplet). 
#Input:LabSystem, theta_values, twotheta, dQ, to_multiplet
function theta_dependence_multiplet(lab::LabSystem,theta_values::Vector{<:Real}, twotheta :: Real, dQ :: Real, to_multiplet::Int64)
    # get multiplets
    energy_values,multiplet_indices=multiplets(lab.eigensys)
    # initalize intensities
    intensities=zeros(length(theta_values))
    # Iterate over all values for dQ
    for i in 1:length(theta_values)
         # set scattering angles
        set_scattering_angles_deg!(lab, theta_values[i],twotheta, dQ)
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