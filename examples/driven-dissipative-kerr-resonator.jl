using LinearAlgebra
using MKL
using DifferentialEquations
using RecursiveArrayTools
using CairoMakie
using Random

include("../matrices.jl")
include("../singlemodeproblemtype.jl")


# define Liouvillian
function liouvillian!(prob::SingleModeCoherentTDVP, α, B, p)
    Δ = p["Δ"]
    κ = p["κ"]
    F = p["F"]
    U = p["U"]
    # Update only necessary operators:
    ad_a!(prob.ad_a, prob.S, α)
    a!(prob.a, prob.S, α)
    adag!(prob.ad, prob.S, α)
    ad2_a2!(prob.ad2_a2, prob.S, α)

    # Your Hamiltonian here:
    prob.H .= Δ*prob.ad_a
    prob.H .+= conj(F)*prob.a
    prob.H .+= F*prob.ad
    prob.H .+= U/2*prob.ad2_a2
    # Update Liouvillian
    lmatrix_efficient!(prob.L, B, prob.S, prob.H, [prob.a], [prob.ad_a], rates=[κ])
end

#time
tmax = 10.0
t_list = range(0, tmax, 1000)

#initial condition for α
α0 = [1.0 + 0.0im]

# System parameters
Δ = 0.1
κ = 1.0
U = 0.5
F = 1.5/sqrt(abs(U))

p = Dict("Δ" => Δ, "κ" => κ, "F" => F, "U" => U)

fig = Figure()
ax = Axis(fig[1, 1])

for order ∈ 1:10
#Initialize TDVP system
TDVPsys = SingleModeCoherentTDVP(order)


#initial conditions
B0 = zeros(ComplexF64, TDVPsys.order, TDVPsys.order)
B0[1,1] = 1
u0 = ArrayPartition(α0, B0)

# Random.seed!(1234)
# α0 = [rand(ComplexF64)]
# B0 = rand(ComplexF64, order, order)
# B0 = B0 + B0'
# u0 = ArrayPartition(α0, B0)

ode_problem = ODEProblem{true}(TDVPsys, u0, (0.0, tmax), p)

sol = solve(ode_problem, Tsit5(), dt = 1e-9, saveat=t_list)
α_vals = [obj.x[1][1] for obj ∈ sol.u]

exp_a = Array{ComplexF64}(undef, length(sol))

for i ∈ eachindex(sol)
    B = sol.u[i].x[2]
    α = sol.u[i].x[1][1]
    smatrix!(TDVPsys.S, α)
    a!(TDVPsys.a, TDVPsys.S, α)
    exp_a[i] = tr(TDVPsys.a[1:end-1, :]*B)/tr(TDVPsys.S[1:end-1, :]*B)
end


# lines!(ax, real.(α_vals), imag.(α_vals))
lines!(ax, real.(exp_a), imag.(exp_a))
end



# Compare to Fock solution
using QuantumOptics

N_fock = 50
basis_Fock = FockBasis(N_fock)
a = destroy(basis_Fock)
ad = create(basis_Fock)
H_Fock = Δ*ad*a + F*(a + ad) + U/2*ad*ad*a*a
ρ0 = dm(coherentstate(basis_Fock, α0[1]))

t_Fock, ρ_Fock = timeevolution.master(t_list, ρ0, H_Fock, [a])

a_Fock = map(x -> expect(a, x), ρ_Fock)

lines!(ax, real.(a_Fock), imag.(a_Fock), linestyle = :dash, color = :black)