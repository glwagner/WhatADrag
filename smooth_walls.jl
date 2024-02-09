using RootSolvers
using GLMakie

struct LaminarLayerLength{FT}
    viscosity :: FT
    coefficient :: FT
end

@inline function (ℓ::LaminarLayerLength)(u★)
    C = ℓ.coefficient
    ν = ℓ.viscosity
    return C * ν / u★
end

function friction_velocity(speed, inner_length_scale;
                           von_karman_constant = 0.4,
                           buoyancy_gradient = 0,
                           reference_height = 10)

    U = speed
    ϰ = von_karman_constant
    z = reference_height
    N = sqrt(buoyancy_gradient)

    function f(u★)
        u★ = max(0, u★)
        ℓ = inner_length_scale(u★)
        h = u★ / N
        z̃ = min(z, h)
        ẑ = max(0, z - h)
        return u★ * log(z̃ / ℓ) + ϰ * N * ẑ - ϰ * U
    end

    method = SecantMethod{Float64}(1e-4, 1e-3)

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

du = 1e-2
U = 1e-1:du:30

viscosity = ν = 1.5e-5
laminar_length_coefficient = Cν = 0.11
charnock_coefficient = α = 0.011
gravitational_acceleration = g = 9.81

fig = Figure()

#axu = Axis(fig[1, 1], xlabel="U(z = 10)", ylabel="u★ / U", xscale=log10)
#axχ = Axis(fig[2, 1], xlabel="U(z = 10)", ylabel="χᴰ = (u★ / U)²", xscale=log10, yscale=log10)
#axℓ = Axis(fig[3, 1], xlabel="U(z = 10)", ylabel="ℓᵤ", xscale=log10, yscale=log10)
axC = Axis(fig[1, 1], xlabel="ū(z = 10 m)", ylabel="Drag coefficient at z = 10 m, Cᴰ(10 m) = (u★ / U)²", xscale=log10, yscale=log10)

ℓν = LaminarLayerLength(ν, Cν)

u★, Cᴰ, ℓ = sweep!(U, ℓν)  
lines!(axC, U, Cᴰ, label="N²=0")

u★, Cᴰ, ℓ = sweep!(U, ℓν, buoyancy_gradient=1e-4)  
lines!(axC, U, Cᴰ, label="N²=10⁻⁴ s⁻²")

u★, Cᴰ, ℓ = sweep!(U, ℓν, buoyancy_gradient=1e-3)  
lines!(axC, U, Cᴰ, label="N²=10⁻³ s⁻²")

axislegend(axC)

#lines!(axu, U, u★ ./ U, label="Neutral, laminar inner layer")
#lines!(axℓ, U, ℓ)

fig

save("stratified_drag_coeff.png", fig)
