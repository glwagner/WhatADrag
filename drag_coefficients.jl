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

function friction_velocity(speed, inner_length_scale;
                           stability_correction = 0,
                           von_karman_constant = 0.4,
                           reference_height = 2)

    U = speed
    κ = von_karman_constant
    h = reference_height

    function f(u★)
        ℓ = inner_length_scale(u★)
        ψ = stability_correction
        return u★ * (log((ℓ + h) / ℓ) - ψ) - κ * U
    end

    method = SecantMethod{Float64}(1e-6, 1e-5)

    sol = find_zero(f, method)

    return sol
end

function sweep!(U, ℓ; kw...)
    N = length(U)
    u★s = zeros(N)
    χs = zeros(N)
    ℓs = zeros(N)
    
    for (n, u) in enumerate(U)
        sol = friction_velocity(u, ℓ; kw...)
        u★ = sol.root
        u★s[n] = u★
        ℓs[n]  = ℓ(u★)
        χs[n]  = (u★ / u)^2
    end

    return u★s, ℓs, χs
end

du = 1e-1
U = du:du:30

viscosity = ν = 1.5e-5
laminar_length_coefficient = Cν = 0.11
charnock_coefficient = α = 0.011
gravitational_acceleration = g = 9.81

fig = Figure()

axu = Axis(fig[1, 1], xlabel="U(z = 2)", ylabel="u★ / U", xscale=log10)
axχ = Axis(fig[2, 1], xlabel="U(z = 2)", ylabel="χᴰ = (u★ / U)²", xscale=log10, yscale=log10)
axℓ = Axis(fig[3, 1], xlabel="U(z = 2)", ylabel="ℓᵤ", xscale=log10, yscale=log10)

ℓν = LaminarLayerLength(ν, Cν)
ℓg = SurfaceWaveLength(g, α)
ℓs = CombinedLengthScale(+, ℓν, ℓg)
ℓm = CombinedLengthScale(max, ℓν, ℓg)

u★, χ, ℓ = sweep!(U, ℓν)  
lines!(axu, U, u★ ./ U, label="Neutral, laminar inner layer")
lines!(axχ, U, χ)
lines!(axℓ, U, ℓ)

u★, χ, ℓ = sweep!(U, ℓg)  
lines!(axu, U, u★ ./ U, label="Neutral, surface wave inner layer")
lines!(axχ, U, χ)
lines!(axℓ, U, ℓ)

u★, χ, ℓ = sweep!(U, ℓs)
lines!(axu, U, u★ ./ U, label="Neutral, laminar + surface wave inner layer")
lines!(axχ, U, χ)
lines!(axℓ, U, ℓ)

#=
u★, χ, ℓ = sweep!(U, ℓm)
lines!(axu, U, u★ ./ U, label="Neutral, max(laminar, surface wave) inner layer")
lines!(axχ, U, χ)
lines!(axℓ, U, ℓ)
=#

u★, χ, ℓ = sweep!(U, ℓs, stability_correction=-5)  
lines!(axu, U, u★ ./ U, label="Stable, surface wave inner layer")
lines!(axχ, U, χ)
lines!(axℓ, U, ℓ)

u★, χ, ℓ = sweep!(U, ℓs, stability_correction=1)  
lines!(axu, U, u★ ./ U, label="Unstable, surface wave inner layer")
lines!(axχ, U, χ)
lines!(axℓ, U, ℓ)

Legend(fig[0, 1], axu)

fig
