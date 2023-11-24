"""
    lmatrix_efficient!(L, B, S, H, [J, JdagJ, Jdag=adjoint.(J), rates=ones(length(J))])

(Efficiently) update the Liouvillian matrix with respect to the Hamiltonian `H`, density matrix `B`, overlap matrix `S` 
and jump operators `J`, `JdagJ`, `Jdag`, with the corresponding rates.

Here all operators are already Kroneckered.
"""
function lmatrix_efficient!(L::T, B::T, S::K, H::T, J::Vector{M}, JdagJ::Vector{M}, orders::Vector{<:Int};
                        Jdag=adjoint.(J), rates=ones(length(J))::Vector{<:Number}) where {T, K<:AbstractMatrix, M<:AbstractMatrix}
    
    size(JdagJ) == size(J) == size(Jdag) == size(rates) || error("J, Jdag, JdagJ, rates do not have the same lengths")
    # size(L) == size(H) == size(S) == (size(B) .+ (1,0)) && all(==(size(H)),size.(J)) == all(==(size(H)),size.(JdagJ)) || error("dimensions of L, H, S, B, J, JdagJ do not match")

    #idcs = CartesianIndices(size(B))
    tmp = Matrix{eltype(L)}(undef, size(B)) # temporary matrix for cache

    # Hamiltonian
    mul!(tmp, B, getindexmat(S,orders), -eltype(B)(im), false) # tmp .= -1im * B * S
    mul!(L, H, tmp, true, false) # L .= H * tmp
    mul!(tmp, B, getindexmat(H,orders), eltype(B)(im), false) # tmp .= 1im * B * H
    mul!(L, S, tmp, true, true) # L .= S * tmp + L

    @inbounds for i ∈ eachindex(J)
        #Quantum Jump: γᵢJᵢ*B*Jᵢ†
        mul!(tmp, B, getindexmatadj(Jdag[i], orders), eltype(B)(rates[i]), false) # tmp .= rates[i] * B * Jdag[i]
        mul!(L, J[i], tmp, true, true) # L .= J[i] * tmp + L

        #Recycling terms: -1/2*γᵢJᵢ†*Jᵢ*B - 1/2*γᵢB*Jᵢ†*Jᵢ
        mul!(tmp, B, getindexmat(JdagJ[i], orders), eltype(B)(-0.5), false) # tmp .= -1/2 * B * JdagJ[i]
        mul!(L, S, tmp, eltype(B)(rates[i]), true) # L .= γᵢ * S * tmp + L

        mul!(tmp, B, getindexmat(S, orders), eltype(B)(-0.5), false) # tmp .= -1/2 * B * S
        mul!(L, JdagJ[i], tmp, eltype(B)(rates[i]), true) # L .= γᵢ * JdagJᵢ * tmp + L
    end
    return L
end




"""
Indexing functions
"""



function getindexmat(S, orders)
    getidx = [trues(o+1, o) for o in orders]
    for mat ∈ getidx
        mat[end,:] .= false
    end
    return reshape(@view(S[kron(getidx...)]), (prod(orders), prod(orders)))
end


# Needed for Jdag
function getindexmatadj(S, orders)
    getidx = [trues(o, o+1) for o in orders]
    for mat ∈ getidx
        mat[:,end] .= false
    end
    return reshape(@view(S[kron(getidx...)]), (prod(orders), prod(orders)))
end



# Now it is correct
# Checked with old implementation
function getindexmat_shiftk(op, orders, k)
    # all getidx are trues except getidx[k] which are falses
    getidx = [trues(o+1, o) for o in orders]
    for i ∈ 1:length(orders)
        i == k ? getidx[i][1,:] .= false : getidx[i][end,:] .= false
    end
   # getidx[k][1:end-1,:] .= false
    #getidx = [trues(o+1, o) for o in orders]
    #getidx[k][end,:] .= true
    #index_mask = [i != k for i in 1:length(orders)]
    return reshape(@view(op[kron(getidx...)]), (prod(orders), prod(orders)))
end






# Tested
function polyalphakron(αk, k, orders)
    # This could perhaps benefit from sparse matrices
    matrices = [Matrix{ComplexF64}(I, n, n) for n ∈ orders]
    matrices[k] .= diagm(orders[k], orders[k], (1=>ones(orders[k]-1)))
    matrices[k][end,:] .= polyα(orders[k], αk)
    return kron(matrices...)
end

# Computes only the column k of the polyalphakron matrix
function polyalphakron_rowk(αk, k, orders)
    # This could perhaps benefit from sparse matrices
    matrices = [Matrix{ComplexF64}(I, n, n) for n ∈ orders]
    matrices[k] = transpose(polyα(orders[k], αk))
    return kron(matrices...)
end


function get_rowk(op, k, orders)
    getidx = [trues(o+1, o) for o in orders]
    for i ∈ 1:length(orders)
        i == k ? getidx[i][1:end-1,:] .= false : getidx[i][end,:] .= false
    end
    index_mask = [i != k for i in 1:length(orders)]
    return reshape(@view(op[kron(getidx...)]), (prod(orders[index_mask]), prod(orders)))
end


function get_row_i_modek(op, k, orders, j)
    getidx = [trues(o+1, o) for o in orders]
    for i ∈ 1:length(orders)
        if i == k
            getidx[i][1:end .!= j,:] .= false
        else 
            getidx[i][end,:] .= false
        end
    end        
    return reshape(@view(op[kron(getidx...)]), (prod(orders[1:end .!= k]), prod(orders)))
end


function get_rowkb(op, k, orders)
    getidx = [trues(o, o) for o in orders]
    for i ∈ 1:length(orders)
        i == k ? getidx[i][1:end-1,:] .= false : nothing
    end
    index_mask = [i != k for i in 1:length(orders)]
    return reshape(@view(op[kron(getidx...)]), (prod(orders[index_mask]), prod(orders)))
end


function get_colkb(op, k, orders)
    getidx = [trues(o, o) for o in orders]
    for i ∈ 1:length(orders)
        i == k ? getidx[i][:,1:end-1] .= false : nothing
    end
    index_mask = [i != k for i in 1:length(orders)]
    return reshape(@view(op[kron(getidx...)]), (prod(orders), prod(orders[index_mask])))
end

function get_col_i_modekb(op, k, orders, j)
    getidx = [trues(o, o) for o in orders]
    for i ∈ 1:length(orders)
        if i == k
            getidx[i][:,1:end .!= j] .= false
        else 
            nothing
        end
    end        
    return reshape(@view(op[kron(getidx...)]), (prod(orders), prod(orders[1:end .!= k])))
end


function get_Skronk(Skron, k, orders)
    getidx = [trues(o+1, o) for o in orders]
    for i ∈ eachindex(orders)
        getidx[i][end,:] .= false
    end
    getidx[k][:,:] .= false
    getidx[k][1,1] = true
    return reshape(@view(Skron[kron(getidx...)]), (prod(orders[1:end .!= k]), prod(orders[1:end .!= k])))
end


function kron_with(S, op, k)
    return reduce(kron, [S[1:k-1]; [op]; S[k+1:end]])
end


function kron_with(S, op1, op2, k, l)
    k < l || error("k must be strictly smaller than l")
    return reduce(kron, [S[1:k-1]; [op1]; S[k+1:l-1]; [op2]; S[l+1:end]])
end