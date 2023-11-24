using LinearAlgebra
using RecursiveArrayTools
using DifferentialEquations
using ProgressLogging
# using Kronecker

include("../matrices.jl")
include("../matricesmultimode.jl")
include("../multimodeproblemtype.jl")

Nmodes = 6
n = 1
orders = n*ones(Int, Nmodes)
prob = CoherentTDVP(orders)

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
        ad2_a2!(prob.a2[k], prob.S[k], α[k])
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
    lmatrix_efficient!(prob.L, B, reduce(kron, prob.S), prob.H, 
                    J, JdagJ, prob.order; rates=κ)
end

# Fscaled = 1.4
F = 1.5695 #sqrt(Nmodes) * Fscaled

# Uscaled = 1
U = 0.1 #Uscaled/Nmodes

Δ = 0.1

J = 0.9

p = Dict("κ" => ones(Nmodes), 
    "F" => F*ones(Nmodes),
    "U" => U*ones(Nmodes),
    "Δ" => Δ*ones(Nmodes),
    "J" => J*ones(Nmodes-1))

tmax = 100.0
t_list = 10 .^ (range(log10(1e-3), log10(tmax), 500))


# Initial state 
α = zeros(ComplexF64, Nmodes) #.+ 1e-12
Bindividual = [zeros(ComplexF64, o, o) for o ∈ prob.order]
for k ∈ eachindex(prob.order)
    Bindividual[k][1,1] = 1
end
B = reduce(kron, Bindividual)

u0 = ArrayPartition(α, B)

# du0 = similar(u0)

# prob(du0, u0, [Δ; κ], nothing)

ode_problem = ODEProblem{true}(prob, u0, (0.0, tmax), p)

sol = solve(ode_problem, Tsit5(), abstol = 1e-6, reltol = 1e-6, dt=1e-8, saveat=t_list, progress = true)




exp_a_3_ord_1 = Array{ComplexF64}(undef, length(sol))

function exp_val(αvec, B, prob, k)
    for j ∈ eachindex(prob.order)
        smatrix!(prob.S[j], αvec[j])
    end
    a!(prob.a[k], prob.S[k], αvec[k])
    return tr(getindexmat(kron_with(prob.S, prob.a[k], k), prob.order)*B) / tr(getindexmat(reduce(kron, prob.S), prob.order) *B)
end

k = 3
for (i, s) ∈ enumerate(sol.u)
    α_i = s.x[1]
    B_i = s.x[2]
    exp_a_3_ord_1[i] = exp_val(α_i, B, prob, k)
end

fig = Figure()
ax = Axis(fig[1,1])
lines!(ax, real.(exp_a_3_ord_1), imag.(exp_a_3_ord_1))
# lines!(ax, real.(exp_a_3_ord_2), imag.(exp_a_3_ord_2))
# scatter!(ax,real(exp_a_3[end])/sqrt(Nmodes), imag.(exp_a_3[end])/sqrt(Nmodes) )