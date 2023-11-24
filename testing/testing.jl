using LinearAlgebra
using MKL
using DifferentialEquations
using RecursiveArrayTools

using CairoMakie

include("../matrices.jl")
include("../singlemodeproblemtype.jl")


# define Liouvillian
function liouvillian!(prob::SingleModeCoherentTDVP, α, B, p)
    Δ = p[1]
    κ = p[2]

    # Update only necessary operators:
    ad_a!(prob.ad_a, prob.S, α)
    a!(prob.a, prob.S, α)

    # Your Hamiltonian here:
    prob.H .= Δ*prob.ad_a
    # Update Liouvillian
    lmatrix_efficient!(prob.L, B, prob.S, prob.H, [prob.a], [prob.ad_a], rates=[κ])
end

order = 20
#Initialize TDVP system
TDVPsys = SingleModeCoherentTDVP(order)

# System parameters
Δ = 2.0
κ = 1.0

#initial conditions
α0 = [1.0 + 0.0im]
B0 = zeros(ComplexF64, TDVPsys.order, TDVPsys.order)
B0[1,1] = 1
u0 = ArrayPartition(α0, B0)

#time
tmax = 10.0
t_list = range(0, tmax, 1000)

ode_problem = ODEProblem{true}(TDVPsys, u0, (0.0, tmax), [Δ, κ])

sol = solve(ode_problem, Tsit5(), abstol = 1e-6, reltol = 1e-6, dt=1e-8, saveat=t_list)
α_vals = [obj.x[1][1] for obj ∈ sol.u]

exp_a = Array{ComplexF64}(undef, length(sol))

for i ∈ eachindex(sol)
    B = sol.u[i].x[2]
    α = sol.u[i].x[1][1]
    smatrix!(TDVPsys.S, α)
    a!(TDVPsys.a, TDVPsys.S, α)
    exp_a[i] = tr(TDVPsys.a[1:end-1, :]*B)/tr(TDVPsys.S[1:end-1, :]*B)
end

fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, real.(α_vals), imag.(α_vals))
lines!(ax, real.(exp_a), imag.(exp_a))