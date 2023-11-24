using LinearAlgebra
using RecursiveArrayTools
using Kronecker

include("../matrices.jl")
include("../matricesmultimode.jl")
include("../multimodeproblemtype.jl")

#Testing
orders = [2,2]
prob = CoherentTDVP(orders)
# α = ComplexF64[1,1]#randn(ComplexF64, 3)

# #randn(ComplexF64, prod(prob.order), prod(prob.order))
# B = zeros(ComplexF64, prod(prob.order), prod(prob.order))
# B[1,1] = 1#randn(ComplexF64, prod(prob.order), prod(prob.order))

# using Random
# Random.seed!(1234)
# # α0 = ComplexF64[1, 1]
# α = randn(ComplexF64, length(orders))
# #B0 = diagm(prod(2 .* ord), prod(2 .* ord),1e-8*ones(prod(2 .* ord)))
# #Bi = diagm(2*ord[1], 2*ord[1], 1e-5*ones(2*ord[1]))
# # Bi = zeros(ComplexF64, ord[1], ord[1])
# # Bi[1,1] =1
# # B0 = kron(Bi, Bi)
# ord = prob.order
# # B01 = randn(ComplexF64, ord[1], ord[1])
# # B01 = B01 + B01'

# # B02 = randn(ComplexF64, ord[2], ord[2])
# # B02 = B02 + B02'
# # B = kron(B01, B02)
# B0 = randn(ComplexF64, prod(ord), prod(ord))
# B= B0 + B0'

function liouvillian!(prob::CoherentTDVP, α, B, p)
    #Decoupled system for checks

    Δ = p[1:2]
    κ = p[3:4]

    # Update only necessary operators:
    for k ∈ eachindex(prob.order)
        ad_a!(prob.ad_a[k], prob.S[k], α[k])
        a!(prob.a[k], prob.S[k], α[k])
    end

    # Your Hamiltonian here:
    # ToDo: Write a helper function to facilitate big operators
    ad_a1 = kronecker(prob.ad_a[1], prob.S[2])
    ad_a2 = kronecker(prob.S[1], prob.ad_a[2])
    a_1 = kronecker(prob.a[1], prob.S[2])
    a_2 = kronecker(prob.S[1], prob.a[2])
    prob.H .= (Δ[1]*ad_a1 .+ Δ[2]*ad_a2)
    # Update Liouvillian
    lmatrix_efficient!(prob.L, B, reduce(kronecker, prob.S), prob.H, 
                    [a_1, a_2], [ad_a1, ad_a2], prob.order;rates=κ)
end




using DifferentialEquations

#time
tmax = 10.0
t_list = range(0, tmax, 1000)

Δ = [2.0, 2.0]
κ = [1.0, 1.0]

# Initial state 
α = ComplexF64[1, 1]
B = zeros(ComplexF64, prod(prob.order), prod(prob.order))
B[1,1] = 1#randn(ComplexF64, prod(prob.order), prod(prob.order))

u0 = ArrayPartition(α, B)

# du0 = similar(u0)

# prob(du0, u0, [Δ; κ], nothing)

ode_problem = ODEProblem{true}(prob, u0, (0.0, tmax), [Δ; κ])

sol = solve(ode_problem, Tsit5(), abstol = 1e-6, reltol = 1e-6, dt=1e-9, saveat=t_list)


# Manual loop
dt = 1e-8
t_current = 0.0

while t_current < tmax
    update!(prob, α)

    # update Liouvillian
    liouvillian!(prob, α, B, p)


    #prob.polα .= transpose(polyα(prob.order, α))

    for k ∈ eachindex(prob.order)
        dα[k] = alpha_derivative(prob, B, k)
    end

    B_derivative!(prob, dB, dα)

    α .+= dα*dt
    B .+= dB*dt

    t_current = t_current + dt
end






















α1_vec = [obj.x[1][1] for obj ∈ sol.u]
α2_vec = [obj.x[1][2] for obj ∈ sol.u]


exp_a_1 = Array{ComplexF64}(undef, length(sol))
exp_a_2 = Array{ComplexF64}(undef, length(sol))


for i ∈ eachindex(sol)
    B = sol.u[i].x[2]
    α1 = sol.u[i].x[1][1]
    α2 = sol.u[i].x[1][2]

    # Update only necessary operators:
    for k ∈ eachindex(prob.order)
        a!(prob.a[k], prob.S[k], sol.u[i].x[1][k])
    end
    a_1 = kron(prob.a[1], prob.S[2])
    a_2 = kron(prob.S[1], prob.a[2])
    S = reduce(kron, prob.S)
    exp_a_1[i] = tr(getindexmat(a_1, prob.order)*B)/tr(getindexmat(S, prob.order)*B)
    exp_a_2[i] = tr(getindexmat(a_2, prob.order)*B)/tr(getindexmat(S, prob.order)*B)
end

fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, real.(α1_vec), imag.(α2_vec))
lines!(ax, real.(exp_a_1), imag.(exp_a_2))
lines!(ax, real.(exp_a_2), imag.(exp_a_2))