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
    
    s1 = sin(θ₁);  s2 = sin(θ₂);  s3 = sin(θ₃)
    s12 = sin(θ₁ - θ₂);  c12 = cos(θ₁ - θ₂)
    s13 = sin(θ₁ - θ₃);  c13 = cos(θ₁ - θ₃)
    s23 = sin(θ₂ - θ₃);  c23 = cos(θ₂ - θ₃)

    A = @SMatrix [
        (m₁+m₂+m₃)*l₁^2    (m₂+m₃)*l₁*l₂*c12    m₃*l₁*l₃*c13;
        (m₂+m₃)*l₁*l₂*c12  (m₂+m₃)*l₂^2         m₃*l₂*l₃*c23;
        m₃*l₁*l₃*c13       m₃*l₂*l₃*c23         m₃*l₃^2
    ]

    b = @SVector [
        (m₂+m₃)*l₁*l₂*s12*ω₂^2 + m₃*l₁*l₃*s13*ω₃^2 - (m₁+m₂+m₃)*g*l₁*s1,
        -(m₂+m₃)*l₁*l₂*s12*ω₁^2 + m₃*l₂*l₃*s23*ω₃^2 - (m₂+m₃)*g*l₂*s2,
        -m₃*l₁*l₃*s13*ω₁^2 - m₃*l₂*l₃*s23*ω₂^2 - m₃*g*l₃*s3
    ]

    return A \ b
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

