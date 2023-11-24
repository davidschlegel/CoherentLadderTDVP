using LinearAlgebra
using RecursiveArrayTools

include("../matrices.jl")
include("../matricesmultimode.jl")
include("../multimodeproblemtype.jl")

#Testing
prob = CoherentTDVP([5, 3, 4])
α = randn(ComplexF64, 3)
ord = prob.order
B01 = randn(ComplexF64, ord[1], ord[1])
B01 = B01 + B01'

B02 = randn(ComplexF64, ord[2], ord[2])
B02 = B02 + B02'
B = kron(B1, B2)
# update S-matrixes
update!(prob, α)
#prob.Skron .= reduce(kron, prob.S)
function liouvillian!(prob, α, B, p)
    prob.L .= randn(ComplexF64, prod(prob.order .+ 1), prod(prob.order))
end


u = ArrayPartition(α, B)
du = similar(u)

alpha_derivative(prob, B, 2)

prob(du, u, nothing, nothing)

# This is still a problemd
# k=1
# Polalphak_row = polyalphakron_rowk(α[k], k, prob.order)

# # Create a view of L of size prod(order)xprod(order)
# L = getindexmat(prob.L, prob.order)

# tmp = Polalphak_row * L

# tmp2 =  prob.polα[k] * [get_row_i_modek(prob.L, k, prob.order, j) for j ∈ 1:prob.order[k]]

# get_row_i_modek(prob.L, k, prob.order, 5)




# Lrowk = get_rowk(prob.L, k, prob.order)



# function efficient_multiply_last_row(A_row, B, orders, k)
#     # Calculate the size of the blocks in B
#     pre_size = prod(orders[1:k-1])
#     post_size = prod(orders[k+1:end])

#     # Reshape B to group the dimensions corresponding to A
#     B_reshaped = reshape(B, (pre_size, orders[k], post_size, pre_size, orders[k], post_size))


#     # Compute only the last row of the result
#     C = zeros(eltype(B), (pre_size^2, orders[k], post_size^2))
#     @inbounds for j ∈ 1:orders[k]
#                 C[:, j, :] .= sum([A_row[i] * B_reshaped[:, i, j, :] for i in 1:orders[k]])
#     end
#     return reshape(C, (prod(orders[1:end .!= k]), prod(orders)))
# end




Polalphak_row = polyalphakron_rowk(α[k], k, prob.order)

# Create a view of L of size prod(order)xprod(order)
L = getindexmat(prob.L, prob.order)





# # # Example values
# m = 3  # Size of A
# N = 4  # Total number of matrices in the Kronecker product
# k = 2  # Position of A in the Kronecker product
# orders = [2, m, 2, 2]  # Sizes of the identity matrices and A in the Kronecker product

# # Example non-trivial matrix A and matrix B
# A = rand(m)  # Non-trivial matrix at position k
# B = rand(prod(orders), prod(orders))  # Matrix to multiply with

# #A = randn(ComplexF64, prob.order[k], prob.order[k])

# lastrow_efficient = efficient_multiply_last_row(A, B, orders, k)

# Amat = zeros(eltype(B), m, m)
# Amat[end,:] .= A
# matrices = [Matrix{Float64}(I, o,o) for o in orders]
# matrices[k] = Amat
# Afull = kron(matrices...)
# productbig = Afull*B