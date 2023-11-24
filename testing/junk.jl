# remove this function
# function polyαb(ord, α)
#     return [binomial(ord, k-1)*α^(ord-k)*(-1)^k for k ∈ 1:ord]
# end


# """
#     lmatrix!(L, B, S, H, J, JdagJ, Jdag=J', rates=ones(length(J)))

# Update the Liouvillian matrix corresponding to the Hamiltonian `H`, density matrix `B`,
# and the jump operators `J`, `JdagJ`, `Jdag`, with the corresponding rates.
# """
# function lmatrix!(L::MArray{Tuple{M, M}, T , 2, K},
#     B::AbstractMatrix{T}, 
#     S::MArray{Tuple{N, N}, T , 2, P},
#     H::MArray{Tuple{N, N}, T , 2, P},
#     J::AbstractVector{MArray{Tuple{N, N}, T , 2, P}},
#     JdagJ::AbstractVector{MArray{Tuple{N, N}, T , 2, P}},
#     Jdag::AbstractVector{MArray{Tuple{N, N}, T , 2, P}} = map(adjoint, J),
#     rates::AbstractVector{<: Real}=ones(length(J)))::Nothing where {N, P, M, K, T <: Complex}

#     size(J) == size(J) == size(Jdag) == size(rates) || error("J, Jdag, JdagJ, rates do not have the same lengths")

#     r = CartesianIndices((N, N) .-1)
#     Heff = H 
#     for (γ,  AdagA) ∈ zip(rates, JdagJ)
#         Heff .-= 1im/2*γ*AdagA
#     end

#     L .= -1im * @view(Heff[r]) * B * @view(S[r])
#     L .+= 1im * @view(S[r]) * B * @view(Heff[r])
#     for (γ, A, Adag) ∈ zip(rates, J, Jdag)
#         L .+= γ * @view(A[r]) * B * @view(Adag[r])
#     end
#     return
# end


# function lmatrix!(L::MArray{Tuple{M, M}, T , 2, K},
#     B::AbstractMatrix{T}, 
#     S::MArray{Tuple{N, N}, T , 2, P},
#     H::MArray{Tuple{N, N}, T , 2, P})::Nothing where {N, P, M, K, T <: Complex}

#     r = CartesianIndices((N, N) .-1)

#     L .= -1im * @view(H[r]) * B * @view(S[r])
#     L .+= 1im * @view(S[r]) * B * @view(H[r])

#     return
# end


# """
# lmatrix_left!(L, B, S, H, [J, JdagJ, Jdag=J', rates=ones(length(J))])

# Update the left derivative of the Liouvillian matrix corresponding to the Hamiltonian `H`, density matrix `B`,
# and the jump operators `J`, `JdagJ`, `Jdag`, with the corresponding rates.
# """
# function lmatrix_left!(L::MArray{Tuple{M, M}, T , 2, K},
#     B::AbstractMatrix{T}, 
#     S::MArray{Tuple{N, N}, T , 2, P},
#     H::MArray{Tuple{N, N}, T , 2, P},
#     J::AbstractVector{MArray{Tuple{N, N}, T , 2, P}},
#     JdagJ::AbstractVector{MArray{Tuple{N, N}, T , 2, P}},
#     Jdag::AbstractVector{MArray{Tuple{N, N}, T , 2, P}} = map(adjoint, J),
#     rates::AbstractVector{<: Real}=ones(length(J)))::Nothing where {N, P, M, K, T <: Complex}

#     size(J) == size(J) == size(Jdag) == size(rates) || error("J, Jdag, JdagJ, rates do not have the same lengths")

#     r = CartesianIndices((N, N) .-1) # indices for no derivative
#     d_l = CartesianIndices((2:N, N-1)) #indices for left derivative
#     Heff = H 
#     for (γ,  AdagA) ∈ zip(rates, JdagJ)
#         Heff .-= 1im/2*γ*AdagA
#     end

#     L .= -1im * @view(Heff[d_l]) * B * @view(S[r])
#     L .+= 1im * @view(S[d_l]) * B * @view(Heff[r])
#     for (γ, A, Adag) ∈ zip(rates, J, Jdag)
#         L .+= γ * @view(A[d_l]) * B * @view(Adag[r])
#     end
#     return
# end


# function lmatrix_left!(L::MArray{Tuple{M, M}, T , 2, K},
#     B::AbstractMatrix{T}, 
#     S::MArray{Tuple{N, N}, T , 2, P},
#     H::MArray{Tuple{N, N}, T , 2, P})::Nothing where {N, P, M, K, T <: Complex}

#     r = CartesianIndices((N, N) .-1) # indices for no derivative
#     d_l = CartesianIndices((2:N, N-1)) #indices for left derivative

#     L .= -1im * @view(H[d_l]) * B * @view(S[r])
#     L .+= 1im * @view(S[d_l]) * B * @view(H[r])

#     return
# end


"""
    smatrix!(S, α)

Update the overlap matrix `S` with coherent field amplitude `α`.
"""
# function smatrix!(Sc::MArray{Tuple{M, N}, T , 2, L}, α::T)::Nothing where {M,N, L, T <: Complex}
#     abs2α = abs2(α)
#     M == N +1 || error("S must be a matrix of size (N+1)x(N)")
#     #We can set the overlap <α|α> = 1, and normalize in the end
#     #Sc[1,1] = 1 # this matrix element never changes
#     for m ∈ 2:M  # border case: first row and column
#         Sc[1,m] = conj(α) * Sc[1,m-1]
#         Sc[m,1] = Sc[1,m]'
#     end
#     Sc[2,2] = abs2α * Sc[1,1] + Sc[1,1]
#     for n ∈ 3:M # border case: second row and column
#         Sc[2,n] = Sc[1,n-1] + abs2α * Sc[1,n-1] + conj(α) * (n-2) * Sc[1,n-2]
#         n == M ? Sc[n,2] = Sc[2,n]' : nothing
#     end
#     for m ∈ 3:M
#         # no conjugation for the diagonal elements
#         Sc[m,m] = Sc[m-1,m-1] + abs2α*Sc[m-1,m-1] + (m-2)*α*Sc[m-2,m-1] + (m-2)*conj(α)*Sc[m-1,m-2] + (m-2)*(m-2)*Sc[m-2,m-2]
#         for n ∈ m+1:M
#             Sc[m,n] = Sc[m-1,n-1] + abs2α*Sc[m-1,n-1] + (m-2)*α*Sc[m-2,n-1] + (n-2)*conj(α)*Sc[m-1,n-2] + (m-2)*(n-2)*Sc[m-2,n-2]
#             # if n == M, do not conjugate
#             n == M ? Sc[n,m] = Sc[m,n]' : nothing 
#         end
#     end
#     return 
# end


# function smatrix!(Sc::MArray{Tuple{2, 1}, T , 2, 2}, α::T)::Nothing where {T <: Complex}
#     Sc[2, 1] = α
#    return
# end


# """
#     smatrix!(α)

# Return the overlap matrix `S` corresponding to coherent field amplitude `α`.
# """
# function smatrix(α::T, ord::N) where {T <: Complex, N<:Integer}
#     # Initialize empty matrix
#     Sc = MMatrix{ord+1, ord, T}(undef)
#     Sc[1,1] = 1
#     # Update matrix
#     smatrix!(Sc, α)
#     return Sc
# end

# function dae_function!(res::AbstractVector{T}, du::AbstractVector{T}, u::AbstractVector{T}, 
#     p::AbstractVector{N}, t::Real) where {T <: Complex, N<:Number}
# Δ = p[1]
# α = u[1]
# B::Matrix{T} = Hermitian(reshape(@view(u[2:end]), (order::Int, order::Int)))
# # update matrices
# smatrix!(S, α)
# ad_a!(ad_a, S, α)

# # Your Hamiltonian here:
# H .= Δ*ad_a

# # update inverse of Sc
# S⁻¹ .= inv(S[r])

# # update Liouvillian
# lmatrix!(L, B, S, H)

# # update left-derivate of Liouvillian
# lmatrix_left!(L_l, B, S, H)

# trY::T = ((transpose(@view(L_l[end, :])) - transpose(@view(S[end, 1:end-1]))*S⁻¹*L) * @view(B[:, end]))[1]
# trC::T = C0_el*(transpose(@view(B[end, :]))*@view(S[r])*@view(B[:, end]))

# # derivative of α
# res[1] = du[1] - (trY + 1e-12) / (trC + 1e-8)

# # drift term
# drift::Matrix{T} = S⁻¹*du[1] * @view(S[r])*B

# # derivative of B
# res[2:end] .= @view(du[2:end]) .- reshape(S⁻¹*L*S⁻¹ - drift - drift', (order::Int)^2)
# end




using Kronecker
using SparseArrays

"""
    c0(S, k)

Returns the matrix C0_k from an array of matrices S representing the overlap matrix in each mode
The resulting array will be of dimension (NxN)^M, where M is the number of modes and N is the dimension of each mode
"""
function c0(S::AbstractArray{<:AbstractMatrix{T}}, k) where T <: Complex
    size_k = size(S[k]) .- 1
    Sk = spzeros(T, size_k)
    Sk[end, end] = factorial(size_k[1])
    Sl = [@view(s[1:end-1, 1:end-1]) for s ∈ S[1:k-1]]
    Sr = [@view(s[1:end-1, 1:end-1]) for s ∈ S[k+1:end]]
    return kronecker(Sl..., Sk, Sr...)
end


"""
    smatrix!(S, α)

Returns the matrices S_i representing the overlap matrix in each mode for α_i 
"""
function smatrix!(S::AbstractVector{<:AbstractMatrix{T}}, α::AbstractVector{T})::Nothing where T <:Complex
    length(S) == length(α) || throw(DimensionMismatch("S and α must have the same length"))
    # Broadcast to update each matrix
    smatrix!.((S, α)...)
    return
end


function kronecker_contract(S, B)
    indices = CartesianIndices(size(B))
    Nmodes = size(S, 1)
    subidcs = CartesianIndices(indices.indices[1:Nmodes])
    C = zeros(eltype(B), size(B))
    # @inbounds for k in eachindex(S)
        @inbounds for i ∈ subidcs
            @inbounds for j ∈ subidcs
                @inbounds for l ∈ subidcs
                    C[CartesianIndex(i, j)] += B[CartesianIndex(i, l)]*prod([S[k][l[k], j[k]] for k ∈ eachindex(S)])
                end
            end
      end
    # end
    return C
end


orders = [3,4,5]
Nmodes = length(orders)
S = [randn(ComplexF64, o, o) for o in orders];
# getidx = [trues(o+1, o) for o in orders]
# for mat ∈ getidx
#    mat[end,:] .= false
# end

# for mat ∈ S
#     mat[end,:] .= 0
# end
B = randn(ComplexF64, (x for y in orders for x in (y, y))...);
#B = randn(ComplexF64, (prod(orders)..., prod(orders)...));
#B = randn(ComplexF64, (orders..., orders...));
B_vec = [randn(ComplexF64, o, o) for o in orders];
B = kronecker(B_vec...);
Bmulti = reshape(B, (orders..., orders...));
Bmulti = permutedims(Bmulti, [1, 4, 2, 5, 3, 6]);
C = kronecker_contract(S, Bmulti);


# Compute the standard Kronecker product
Skron = kron(S...)
# Reshape B appropriately
# B_reshaped = reshape(permutedims(B, [1, 4, 2, 5, 3, 6]), prod(orders), prod(orders))
#permutedims(reshape(B, prod(orders), prod(orders)), (2, 1))
#B_reshaped = reshape(B, prod(orders), prod(orders))
# Perform matrix multiplication
C_standard = B*Skron

# Reshape the result to compare
C_standard_reshaped = reshape(C_standard, size(B))

# Compare the results
difference = norm(C - C_standard_reshaped)
println("Difference between the methods: $difference")




Skron = kron(S...);
B_big = reshape(B, (prod(orders), prod(orders)));
C2 = B_big*Skron;

C2 ≈ reshape(C, (prod(orders), prod(orders)))


    # @inbounds for k ∈ eachindex(prob.order)
    #     tmp = zeros(T, size(prob.drift))
    #     @inbounds for col ∈ 1:prob.order[k]
    #         tmp .+= get_col_i_modekb(B, k, prob.order, col) * prob.polα[k][col]
    #     end
    #     # only add to all columns k (except first column)
    #     get_colk_shiftedb(tmp, k, prob.order) .+= get_colkb(B, k, prob.order)
    #     #multiply everything by conj(dα[k])
    #     tmp .*= conj(dα[k])

    #     # add tmp to drift
    #     if k == 1
    #         prob.drift .= tmp
    #     else
    #         prob.drift .+= tmp
    #     end
    # end
    #     # B*polyα_k
    #     tmp .+= prob.polα[k][j] * get_row_i_modek(prob.L, k, prob.order, j)
    # conj(dα[k])

    # # drift term, avoids inverse
    # prob.drift .= B[:,end] * prob.polα 
    # prob.drift[:, 2:end] .+= B[:, 1:end-1]
    # prob.drift .*= conj(du.x[1][1])


    # function get_all_exceptk(op, k, orders)
#     getidx = [Matrix(I, o,o) for o in orders]
#     getidx[k][:,:] .= true
#     return reshape(@view(op[kron(getidx...)]), (orders[k], orders[k]))
# end


#Testing lmatrix_efficient! To be removed later
# orders = [3,4]
# S = [randn(ComplexF64, o+1, o) for o in orders];
# B = randn(ComplexF64, prod(orders), prod(orders));
# Skron = reduce(kron, S)
# H = randn(ComplexF64, size(Skron))
# L = similar(H)

# J =  [randn(ComplexF64, size(Skron))]
# JdagJ = [randn(ComplexF64, size(Skron))]
# rates = [rand(ComplexF64)]

# lmatrix_efficient!(L, B, Skron, H, J, JdagJ, orders; rates=rates)

# k=2
# Lbk = getindexmat_shiftk(L, orders, k)
# Sinv = kron([inv(S[i][1:end-1,:]) for i in 1:length(orders)]...)
# S_l = getindexmat_shiftk(Skron, orders, k)

# #Test setup
# orders = [3,4,15]
# α = [randn(ComplexF64) for i ∈ 1:length(orders)]

# S = [randn(ComplexF64, i+1, i) for i ∈ orders]
# for i ∈ eachindex(orders)
#     S[i][1,1] = 1
#     smatrix!(S[i], α[i])
# end

# Skron = kron(S...)
# Sinv = kron([inv(Hermitian(S[i][1:end-1,:])) for i in 1:length(orders)]...)
# S_l = getindexmat_shiftk(Skron, orders, k)
# S_polk = polyalphakron(α[k], k, orders)
# H = randn(ComplexF64, size(Skron))
# L = similar(H)
# J = [ randn(ComplexF64, size(Skron)) for i ∈ 1:2]
# JdagJ = [ randn(ComplexF64, size(Skron)) for i ∈ 1:2]
# rates = [2.0, 1.5]
# B = randn(ComplexF64, prod(orders), prod(orders))
# lmatrix_efficient!(L, B, Skron, H, J, JdagJ, orders; rates=rates)


# #Y0k
# k=2
# Lbk = getindexmat_shiftk(L, orders, k)
# Spolα = polyalphakron(α[k], k, orders)
# Y0k = Lbk - Spolα*getindexmat(L, orders)

# #Get column of Y0k
# Lbk_row = get_rowk(L, k, orders)
# S_lSinv_L =  Spolα*getindexmat(L, orders)
# S_lSinv_L_row = get_rowkb(S_lSinv_L, k, orders)
# Y0k_row = Lbk_row - S_lSinv_L_row


# Bk_col = get_colkb(B, k, orders)


# tr_easy = tr(Y0k_row*Bk_col)
# tr_hard = tr(Y0k*B)


# Polalphak_row = polyalphakron_rowk(α[k], k, orders)
# Polalphak = polyalphakron(α[k], k, orders)
# Polalphak_row ≈ get_rowkb(Polalphak, k, orders)



# S_lSinv_L_rowb = Polalphak_row*getindexmat(L, orders)
# S_lSinv_L_row ≈ S_lSinv_L_colb


# Y_col_easy = Lbk_row - S_lSinv_L_rowb


# #C0
# c0el = factorial(orders[k])
# Skronk = reduce(kron, [mat[1:end-1, :] for mat ∈ S[1:end .!= k]])
# Skronk_from_Skron = get_Skronk(Skron, k, orders)

# Skronk ≈ Skronk_from_Skron
# Bk_row = get_rowkb(B, k, orders)
# Bk_col = get_colkb(B, k, orders)

# C = c0el*Skronk_from_Skron*Bk_row*getindexmat(Skron, orders)*Bk_col

# trC = tr(C)



# # Calculates an indexing function to index Kroneckered matrices
# # k = 0: take the non-shifted Block, i.e. the upper left corner
# # k > 0: shift Block of mode k down.
# #       Corresponds to the left derivative.
# # revision: done
# function index(ord, k)
#     matrices = Array{BitMatrix}(undef, length(ord))
#     for (j,i) ∈ enumerate(ord)
#         x = falses(i+1, i)
#         j==k ? x[2:end, :] = trues(i, i) : x[1:end-1, :] = trues(i, i)
#         matrices[j] = x
#     end
#     return kron(matrices...)
# end


# idcsk = [index(orders,k) for k in 1:length(orders)] #Kroneckered Bit-Matrix where mode k is block shifted.
# idl(x, k) = reshape(x[idcsk[k]], prod(orders), prod(orders))
# S_l_old = idl(Skron, k)