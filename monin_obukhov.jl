using RootSolvers
using GLMakie

struct LaminarLayerLength{FT}
    viscosity :: FT
    coefficient :: FT
end

struct SurfaceWaveLength{FT}
    gravitational_acceleration :: FT
    coefficient :: FT
end

@inline function (ℓ::LaminarLayerLength)(u★)
    C = ℓ.coefficient
    ν = ℓ.viscosity
    return C * ν / u★
end

@inline function (ℓ::SurfaceWaveLength)(u★) 
    α = ℓ.coefficient
    g = ℓ.gravitational_acceleration
    return α * u★^2 / g
end

struct CombinedLengthScale{N, O, L}
    op :: O
    length_scales :: L

    function CombinedLengthScale(op, length_scales...)
        length_scales = tuple(length_scales...)
        N = length(length_scales)
        O = typeof(op)
        L = typeof(length_scales)
        return new{N, O, L}(op, length_scales)
    end
end

@inline function (cls::CombinedLengthScale{N})(u★) where N
    ℓs = cls.length_scales    
    ℓ₁ = ℓs[1]
    s = ℓ₁(u★)
    ntuple(Val(N-1)) do n
        ℓ⁺ = ℓs[n+1](u★)
        s = cls.op(s, ℓ⁺)
    end
    return s
end

struct MoistAirState{Q, T, U, V}
    q :: Q
    θ :: T
    u :: U
    v :: V
end

MoistAirState(q, θ, u) = MoistAirState(q, θ, u, nothing)

Base.@kwdef struct SimilarityParameters{FT}
    gravitational_acceleration :: FT = 9.81
    von_karman_constant :: FT = 0.4
    molar_mass_ratio = 1.61
end

struct UniversalFunction{FT}
    a :: FT
    b :: FT
    c :: FT
end

@inline function (ψ::UniversalFunction)(Ri)
    a = ψ.a
    b = ψ.b
    c = ψ.c

    ϕ⁻¹ = (1 - b * Ri)^c
    ψ_unstable = log((1 + ϕ⁻¹)^2 * (1 + ϕ⁻¹^2) / 8) - 2 * atan(ϕ⁻¹) + π/2
    ψ_stable = - a * Ri
    return ifelse(Ri < 0, ψ_unstable, ψ_stable)
end



#=
Base.@kwdef struct Businger{FT}
    au :: FT = 4.7
    bu :: FT = 15.0
    cu :: FT = 0.25
    ac :: FT = 6.4
    bc :: FT = 9.0
    cc :: FT = 0.5
end
=#

@inline function ψ(a, b, c, Ri)
    ϕ⁻¹ = (1 - b * Ri)^c
    ψ_unstable = log((1 + ϕ⁻¹)^2 * (1 + ϕ⁻¹^2) / 8) - 2 * atan(ϕ⁻¹) + π/2

    ψ_stable = - a * Ri

    return ifelse(Ri < 0, ψ_unstable, ψ_stable)
end


@inline function ψu(func::Businger, Ri)
    a = func.au
    b = func.bu
    c = func.cu

    ϕ⁻¹          = (1 - b * Ri)^c
    ψ_convecting = log((1 + ϕ⁻¹)^2 * (1 + ϕ⁻¹^2) / 8) - 2 * atan(ϕ⁻¹) + π/2
    ψ_stable     = - a * Ri

    return ifelse(Ri < 0, ψ_convecting, ψ_stable)
end

@inline function ψc(func::Businger, Ri)
    a = func.ac
    b = func.bc
    c = func.cc

    ϕ⁻¹          = (1 - b * Ri)^c
    ψ_convecting = 2 * log((1 + ϕ⁻¹) / 2)
    ψ_stable     = - a * Ri

    return ifelse(Ri < 0, ψ_convecting, ψ_stable)
end
    
@inline momentum_flux_scale(f, h, ℓ, Ri) = 1 / (log(h/ℓ) - ψu(f, Ri) + ψu(f, ℓ * Ri / h))
@inline   tracer_flux_scale(f, h, ℓ, Ri) = 1 / (log(h/ℓ) - ψc(f, Ri) + ψc(f, ℓ * Ri / h))

@inline function virtual_temperature(parameters, state)
    r = parameters.molar_mass_ratio
    δ = r - 1
    q = state.q
    θ = state.θ
    return θ * (1 + δ * q)
end

function buoyancy_scale(θ★, q★, surface_state, parameters)
    θ★ = fluxes.θ
    q★ = fluxes.q
    𝒯₀ = virtual_temperature(parameters, surface_state)
    q₀ = surface_state.q
    θ₀ = surface_state.θ
    r = parameters.molar_mass_ratio
    g = parameters.gravitational_acceleration
    δ = r - 1
    b★ = g / 𝒯₀ * (θ★ * (1 + δ * q₀) + δ * θ₀ * q★)
    return b★
end

function fixed_point_fluxes(u★, θ★, q★,
                            surface_state,
                            inner_length_scales,
                            universal_function,
                            parameters) 

    Δu = differences.u
    Δv = differences.v
    Δθ = differences.θ
    Δq = differences.q

    ϰ = parameters.von_karman_constant
    f = universal_function

    b★ = buoyancy_scale(θ★, q★, surface_state, parameters)
    Riₕ = - ϰ * h * b★ / u★^2

    ℓu = inner_length_scales.u(u★)
    ℓθ = inner_length_scales.θ(u★)
    ℓq = inner_length_scales.q(u★)

    χu = momentum_flux_scale(f, h, ℓu, Riₕ)
    χθ =   tracer_flux_scale(f, h, ℓθ, Riₕ)
    χq =   tracer_flux_scale(f, h, ℓq, Riₕ)

    u★ = ϰ * χu * sqrt(Δu^2 + Δv^2)
    θ★ = ϰ * χθ * Δθ
    q★ = ϰ * χq * Δq

    return u★, θ★, q★
end

@inline function fluxes_not_converged(δu::FT, δθ, δq) where FT
    r = 1e-5
    return (δu > r) | (δθ > r) | (δq > r)
end

function compute_fluxes(::Val{N},
                        fluxes,
                        differences,
                        surface_state,
                        inner_length_scales,
                        reference_height,
                        universal_function,
                        parameters) where N

    u★ = fluxes.u
    θ★ = fluxes.θ
    q★ = fluxes.q

    δu = Inf
    δθ = Inf
    δq = Inf
    
    iter = 1

    while fluxes_not_converged(δu, δθ, δq)
        u★⁺, θ★⁺, q★⁺ = fixed_point_fluxes(u★, θ★, q★,
                                           surface_state,
                                           inner_length_scales,
                                           universal_function,
                                           parameters) 

        δu = u★⁺ - u★
        δθ = θ★⁺ - θ★
        δq = q★⁺ - q★

        u★ = u★⁺
        θ★ = θ★⁺
        q★ = q★⁺ 

        @show iter
        @show u★ θ★ q★
        @show δu δθ δq

        iter += 1
    end                         
end

# Test
surface_state = MoistAirState(0.02, 20, 0, 0)
differences   = MoistAirState(0.01,  5, -10, 20)
parameters    = SimilarityParameters() 
h = 2
g = parameters.gravitational_acceleration
α = 0.011

ℓu = SurfaceWaveLength(g, α)
ℓθ = SurfaceWaveLength(g, α)
ℓq = SurfaceWaveLength(g, α)

length_scales = MoistAirState(ℓq, ℓθ, ℓu)

func = Businger()

fluxes = MoistAirState(1e-2, 1e-2, 1e-2)

sol = compute_fluxes(Val(20), fluxes, differences, surface_state,
                     length_scales, h, func, parameters)

