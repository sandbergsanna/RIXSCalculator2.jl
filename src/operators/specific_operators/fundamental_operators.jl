################################################################################
#
#   INDIVIDUAL OPERATORS FOR L AND S
#   PURELY BASED ON THEIR NATURAL QUANTUM NUMBERS
#
################################################################################

# delta function
"""
    delta(i,j)

Standard Kroenecker delta ``\\delta_{i,j}``.
"""
function delta(i, j)
    if i == j
        return 1
    else
        return 0
    end
end

# dagger of something
"""
    dagger(x)

The function computes the complex conjugate of a number, ``x^\\dagger``.
"""
function dagger(x)
    return x'
end


# export delta and dagger
export delta, dagger





# operator functions for L
# corresponding to matrix elements
# < l,mlp | operator | l,ml >
"""
    operatorLz(l::Int64, mlp::Int64, ml::Int64) :: Complex{Float64}

This function corresponds to the matrix element ``\\left< l,m_{lp} \\left| L_z \\right| l,m_l \\right>.
"""
function operatorLz(l::Int64, mlp::Int64, ml::Int64) :: Complex{Float64}
    return delta(mlp, ml)   * ml
end

"""
    operatorLplus(l::Int64, mlp::Int64, ml::Int64) :: Complex{Float64}

This function corresponds to the matrix element ``\\left< l,m_{lp} \\left| L^+ \\right| l,m_l \\right>.
"""
function operatorLplus(l::Int64, mlp::Int64, ml::Int64) :: Complex{Float64}
    return delta(mlp, ml+1) * sqrt(l*(l+1) - ml*(ml+1))
end

"""
    operatorLminus(l::Int64, mlp::Int64, ml::Int64) :: Complex{Float64}

This function corresponds to the matrix element ``\\left< l,m_{lp} \\left| L^- \\right| l,m_l \\right>.
"""
function operatorLminus(l::Int64, mlp::Int64, ml::Int64) :: Complex{Float64}
    return delta(mlp, ml-1) * sqrt(l*(l+1) - ml*(ml-1))
end

"""
    operatorLx(l::Int64, mlp::Int64, ml::Int64) :: Complex{Float64}

This function corresponds to the matrix element ``\\left< l,m_{lp} \\left| L_x \\right| l,m_l \\right>.
"""
function operatorLx(l::Int64, mlp::Int64, ml::Int64) :: Complex{Float64}
    return 0.5 * (operatorLplus(l,mlp,ml) + operatorLminus(l,mlp,ml))
end

"""
    operatorLy(l::Int64, mlp::Int64, ml::Int64) :: Complex{Float64}

This function corresponds to the matrix element ``\\left< l,m_{lp} \\left| L_y \\right| l,m_l \\right>.
"""
function operatorLy(l::Int64, mlp::Int64, ml::Int64) :: Complex{Float64}
    return -0.5*im * (operatorLplus(l,mlp,ml) - operatorLminus(l,mlp,ml))
end



# export operators
export operatorLx, operatorLy, operatorLz, operatorLplus, operatorLminus


# operator functions for S
# corresponding to matrix elements
# < s,msp | operator | s,ms >
"""
    operatorSz(s::Rational, msp::Rational, ms::Rational) :: Complex{Float64}

This function corresponds to the matrix element ``\\left< s,m_{sp} \\left| S_z \\right| s,m_s \\right>.
"""
function operatorSz(s::Rational, msp::Rational, ms::Rational) :: Complex{Float64}
    return delta(msp, ms)   * ms
end

"""
    operatorSplus(s::Rational, msp::Rational, ms::Rational) :: Complex{Float64}

This function corresponds to the matrix element ``\\left< s,m_{sp} \\left| S^+ \\right| s,m_s \\right>.
"""
function operatorSplus(s::Rational, msp::Rational, ms::Rational) :: Complex{Float64}
    return delta(msp, ms+1) * sqrt(s*(s+1) - ms*(ms+1))
end

"""
    operatorSminus(s::Rational, msp::Rational, ms::Rational) :: Complex{Float64}

This function corresponds to the matrix element ``\\left< s,m_{sp} \\left| S^- \\right| s,m_s \\right>.
"""
function operatorSminus(s::Rational, msp::Rational, ms::Rational) :: Complex{Float64}
    return delta(msp, ms-1) * sqrt(s*(s+1) - ms*(ms-1))
end

"""
    operatorSx(s::Rational, msp::Rational, ms::Rational) :: Complex{Float64}

This function corresponds to the matrix element ``\\left< s,m_{sp} \\left| S_x \\right| s,m_s \\right>.
"""
function operatorSx(s::Rational, msp::Rational, ms::Rational) :: Complex{Float64}
    return 0.5 * (operatorSplus(s,msp,ms) + operatorSminus(s,msp,ms))
end

"""
    operatorSy(s::Rational, msp::Rational, ms::Rational) :: Complex{Float64}

This function corresponds to the matrix element ``\\left< s,m_{sp} \\left| S_y \\right| s,m_s \\right>.
"""
function operatorSy(s::Rational, msp::Rational, ms::Rational) :: Complex{Float64}
    return -0.5*im * (operatorSplus(s,msp,ms) - operatorSminus(s,msp,ms))
end

# export operators
export operatorSx, operatorSy, operatorSz, operatorSplus, operatorSminus