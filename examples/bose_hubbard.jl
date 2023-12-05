using LinearAlgebra
using RecursiveArrayTools
using DifferentialEquations
using ProgressLogging
using CairoMakie
# using Kronecker

include("../matrices.jl")
include("../matricesmultimode.jl")
include("../multimodeproblemtype.jl")



function liouvillian!(prob::CoherentTDVP, α, B, p)
    #Decoupled system for checks

    J = p["J"]
    κ = p["κ"]
    F = p["F"]
    U = p["U"]
    Δ = p["Δ"]

    # Update only necessary operators:
    for k ∈ eachindex(prob.order)
        a!(prob.a[k], prob.S[k], α[k])
        adag!(prob.ad[k], prob.S[k], α[k])
        ad_a!(prob.ad_a[k], prob.S[k], α[k])
        ad2_a2!(prob.ad2_a2[k], prob.S[k], α[k])
    end

    prob.H .= 0

    for k ∈ eachindex(prob.order)
                prob.H .-= Δ[k]*kron_with(prob.S, prob.ad_a[k], k)
                prob.H .+= U[k]/2*kron_with(prob.S, prob.ad2_a2[k], k)
                prob.H .+= F[k]*kron_with(prob.S, prob.a[k], k)
                prob.H .+= conj(F[k])*kron_with(prob.S, prob.ad[k], k)
                # loop over nearest neighbors
                if k < length(prob.order)
                    prob.H .-= J[k]*kron_with(prob.S, prob.a[k], prob.ad[k+1], k, k+1)
                    prob.H .-= J[k]*kron_with(prob.S, prob.ad[k], prob.a[k+1], k, k+1)
                end
    end
    J = [kron_with(prob.S, prob.a[k], k) for k ∈ eachindex(prob.order)]
    JdagJ = [kron_with(prob.S, prob.ad_a[k], k) for k ∈ eachindex(prob.order)]
    lmatrix_efficient!(prob.L, B, prob.Skron, prob.H, 
                    J, JdagJ, prob.order; rates=κ)
end


Nmodes = 3
n = 2
orders = n*ones(Int, Nmodes)
prob = CoherentTDVP(orders)

Fscaled = 2.2
# F = sqrt(Nmodes) * Fscaled

Uscaled = 1.0
# U = 0.1 #Uscaled/Nmodes

Δ = 3.0

J = 1.0*2

p = Dict("κ" => ones(Nmodes), 
    "F" => Fscaled*ones(Nmodes),
    "U" => Uscaled*ones(Nmodes),
    "Δ" => Δ*ones(Nmodes),
    "J" => J*ones(Nmodes-1))

tmax = 100.0
t_list = 10 .^ (range(log10(1e-3), log10(tmax), 500))


# Initial state 
# α = zeros(ComplexF64, Nmodes) #.+ 1e-12
α = sol.u[end].x[1]
Bindividual = [zeros(ComplexF64, o, o) for o ∈ prob.order]
for k ∈ eachindex(prob.order)
    Bindividual[k][1,1] = 1
end

B = reduce(kron, Bindividual)

u0 = ArrayPartition(α, B)

# du0 = similar(u0)

# prob(du0, u0, [Δ; κ], nothing)

ode_problem = ODEProblem{true}(prob, u0, (0.0, tmax), p)

sol = solve(ode_problem, Tsit5(), abstol = 1e-6, reltol = 1e-6, dt=1e-8, saveat=t_list, progress = true, callback = cb)

# Define a condition for the callback (e.g., at every time step)
condition(u, t, integrator) = true  # This condition is always true; adjust as needed

# Define the action to be taken when the condition is met
function affect!(integrator)
    println("Time: ", integrator.t, " | α: ", integrator.u.x[1])
    # Additional actions can be added here
end

# Create the callback
cb = DiscreteCallback(condition, affect!)

exp_a_3_ord_1 = Array{ComplexF64}(undef, length(sol))

function exp_val(αvec, B, prob, k)
    for j ∈ eachindex(prob.order)
        smatrix!(prob.S[j], αvec[j])
    end
    a!(prob.a[k], prob.S[k], αvec[k])
    return tr(getindexmat(kron_with(prob.S, prob.a[k], k), prob.order)*B) / tr(getindexmat(reduce(kron, prob.S), prob.order) *B)
end

k = 1
for (i, s) ∈ enumerate(sol.u)
    α_i = s.x[1]
    B_i = s.x[2]
    exp_a_3_ord_1[i] = exp_val(α_i, B, prob, k)
end

using LaTeXStrings
fig = Figure()
ax = Axis(fig[1,1],
    xlabel = L"\mathrm{Re}(\langle \hat{a}_3\rangle)",
    ylabel = L"\mathrm{Im}(\langle \hat{a}_3\rangle)",
    title = "Bose-Hubbard model",
    backgroundcolor = :transparent)
lines!(ax, real.(exp_a_3_ord_1), imag.(exp_a_3_ord_1), label = L"n=2")
# lines!(ax, real.(exp_a_3_ord_2), imag.(exp_a_3_ord_2), label = L"n=2")

#make Legend
axislegend(ax)

#Show information in a box
# Number of modes, F, U, Δ, J
text!(ax, 0.0, 1.0; 
    text="N = $Nmodes \nF = $Fscaled \nU = $Uscaled \nΔ = $Δ \nJ = $J",
    align = (:left, :top))
# N = $Nmodes \n F = $F \n U = $U \n Δ = $Δ \n J = $J