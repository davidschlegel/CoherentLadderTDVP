include("matrices.jl")
include("odefunction.jl")
using DifferentialEquations
using CairoMakie

# pre-allocate memory for the TDVP algorithm
order::Int = 10
S = smatrix(1.0+0.0im, order)
ad_a = MMatrix{order+1, order, ComplexF64}(undef)
H = MMatrix{order+1, order, ComplexF64}(undef)
S⁻¹ = MMatrix{order, order, ComplexF64}(undef)
L = MMatrix{order+1, order, ComplexF64}(undef)

C0_el = factorial(order)


# system parameters
const Δ = 2.0
p = [Δ]

# initial conditions
α0 = 1.0 + 0.0im
B0 = zeros(ComplexF64, order^2)
B0[1] = 1
u0 = [α0; B0]

du = zeros(ComplexF64, order^2+1)#Array{ComplexF64}(undef, order^2+1)

# call once to compile
ode_function!(du, u0, [Δ], 0.0)


# time 
tmax = 10
nstep = 1000
tspan = (0, tmax)
t_list = range(0, tmax, nstep)

# ODE Problem
prob = ODEProblem{true}(ode_function!, u0, tspan, p, saveat=t_list)
sol = solve(prob, Tsit5(), dt=1e-8)
# abstol=1e-6, reltol=1e-6

#prob = DAEProblem{true}(dae_function!, du, u0, tspan, p, saveat=t_list, differential_vars=trues(size(u0)))
#sol = solve(prob, DFBDF(autodiff=false), initializealg=NoInit(), dt = 1e-8)#initializealg=BrownFullBasicInit()

#epxectation value
exp_a = Array{ComplexF64}(undef, length(sol))
a = MMatrix{order+1, order, ComplexF64}(undef)

for i ∈ eachindex(sol)
    B = reshape(sol[2:end, i], order, order)
    α = sol[1, i]
    smatrix!(S, α)
    a!(a, S, α)
    exp_a[i] = tr(a[1:end-1, :]*B)/tr(S[1:end-1, :]*B)
end