using StaticArrays

mutable struct Pendulum{N,T<:Real}
    M :: SVector{N,T}
    L :: SVector{N,T}
    θ :: MVector{N,T}
    ω :: MVector{N,T}
end

function Pendulum(n; M::AbstractVector{T}=ones(n), L::AbstractVector{T}=ones(n), θ::AbstractVector{T}=zeros(n), ω::AbstractVector{T}=zeros(n)) where {T<:Real}
    @assert all(length(x)==n for x in (M,L,θ,ω)) "length mismatch"
    Pendulum{n,T}(SVector{n,T}(M), SVector{n,T}(L), MVector{n,T}(θ), MVector{n,T}(ω))
end

Pendulum(M, L, θ, ω) = Pendulum(n=length(M); M=M, L=L, θ=θ, ω=ω)

function RK4(state, h, derivatives_func, params...)
    k1 = h * derivatives_func(state, params...)
    k2 = h * derivatives_func(state + k1/2, params...)
    k3 = h * derivatives_func(state + k2/2, params...)
    k4 = h * derivatives_func(state + k3, params...)
    return state + (k1 + 2k2 + 2k3 + k4) / 6
end

function pendulum_derivatives(state, acceleration_func, params...)
    θs, ωs = state
    dθs = ωs
    dωs = acceleration_func(θs, ωs, params...)
    return [dθs, dωs]
end

function single_pendulum_accelerations(θs, ωs, masses, lengths, g)
    θ₁ = θs[1]
    l₁ = lengths[1]
    
    dω₁ = -(g/l₁) * sin(θ₁)
    return [dω₁]
end

function double_pendulum_accelerations(θs, ωs, masses, lengths, g)
    θ₁, θ₂ = θs
    ω₁, ω₂ = ωs
    m₁, m₂ = masses
    l₁, l₂ = lengths

    Δθ = θ₁ - θ₂
    sin_Δθ, cos_Δθ = sincos(Δθ)
    
    dω₁ = (-g * (2*m₁ + m₂) * sin(θ₁) -
           m₂ * g * sin(θ₁ - 2*θ₂) -
           2 * sin_Δθ * m₂ * (ω₂^2 * l₂ + ω₁^2 * l₁ * cos_Δθ)) /
          (l₁ * (2*m₁ + m₂ - m₂ * cos(2*Δθ)))
    
    dω₂ = (2 * sin_Δθ * (ω₁^2 * l₁ * (m₁ + m₂) +
           g * (m₁ + m₂) * cos(θ₁) +
           ω₂^2 * l₂ * m₂ * cos_Δθ)) /
          (l₂ * (2*m₁ + m₂ - m₂ * cos(2*Δθ)))
    return [dω₁, dω₂]
end

function triple_pendulum_accelerations(θs, ωs, masses, lengths, g)
    θ₁, θ₂, θ₃ = θs
    ω₁, ω₂, ω₃ = ωs
    m₁, m₂, m₃ = masses
    l₁, l₂, l₃ = lengths
    
    # Precompute trig functions
    s₁, c₁ = sincos(θ₁)
    s₂, c₂ = sincos(θ₂)  
    s₃, c₃ = sincos(θ₃)
    s₁₂, c₁₂ = sincos(θ₁ - θ₂)
    s₁₃, c₁₃ = sincos(θ₁ - θ₃)
    s₂₃, c₂₃ = sincos(θ₂ - θ₃)
    
    # Mass matrix (inertia matrix)
    M11 = (m₁ + m₂ + m₃) * l₁^2
    M12 = (m₂ + m₃) * l₁ * l₂ * c₁₂
    M13 = m₃ * l₁ * l₃ * c₁₃
    M21 = (m₂ + m₃) * l₁ * l₂ * c₁₂
    M22 = (m₂ + m₃) * l₂^2
    M23 = m₃ * l₂ * l₃ * c₂₃
    M31 = m₃ * l₁ * l₃ * c₁₃
    M32 = m₃ * l₂ * l₃ * c₂₃
    M33 = m₃ * l₃^2
    
    # Coriolis and centrifugal terms
    h₁ = -(m₂ + m₃) * l₁ * l₂ * s₁₂ * ω₂^2 - m₃ * l₁ * l₃ * s₁₃ * ω₃^2
    h₂ = (m₂ + m₃) * l₁ * l₂ * s₁₂ * ω₁^2 - m₃ * l₂ * l₃ * s₂₃ * ω₃^2
    h₃ = m₃ * l₁ * l₃ * s₁₃ * ω₁^2 + m₃ * l₂ * l₃ * s₂₃ * ω₂^2
    
    # Gravity terms
    G₁ = -(m₁ + m₂ + m₃) * g * l₁ * s₁
    G₂ = -(m₂ + m₃) * g * l₂ * s₂
    G₃ = -m₃ * g * l₃ * s₃
    
    # Right hand side vector
    rhs₁ = h₁ + G₁
    rhs₂ = h₂ + G₂ 
    rhs₃ = h₃ + G₃
    
    # Solve M * accelerations = rhs
    M_matrix = @SMatrix [
        M11 M12 M13;
        M21 M22 M23;
        M31 M32 M33
    ]
    
    rhs = @SVector [rhs₁, rhs₂, rhs₃]
    
    return M_matrix \ rhs
end


function (pen::Pendulum{1})(; dt=0.01, g=9.81)
    pen.θ, pen.ω = RK4([pen.θ, pen.ω], dt, pendulum_derivatives, single_pendulum_accelerations, pen.M, pen.L, g)
    return pen
end

function (pen::Pendulum{2})(; dt=0.01, g=9.81)
    pen.θ, pen.ω = RK4([pen.θ, pen.ω], dt, pendulum_derivatives, double_pendulum_accelerations, pen.M, pen.L, g)
    return pen
end

function (pen::Pendulum{3})(; dt=0.01, g=9.81)
    pen.θ, pen.ω = RK4([pen.θ, pen.ω], dt, pendulum_derivatives, triple_pendulum_accelerations, pen.M, pen.L, g)
    return pen
end

