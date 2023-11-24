const t = transpose


"""
    smatrix!(S, α)

Update the overlap matrix `S` with coherent field amplitude `α`.
Note that the matrix S is of dimension (N+1)xN.
"""
function smatrix!(Sc::AbstractMatrix{T}, α::T)::Nothing where {T <: Complex}
    abs2α = abs2(α)
    size(Sc)[1] - 1 == size(Sc)[2] || error("S must be a matrix of size (N+1)x(N)")
    # We can set the overlap <α|α> = 1, and normalize in the end
    # Sc[1,1] = 1 # this matrix element never changes, so it is better to have it initialized already
    
    # special case 2 x 1 matrix
    M = size(Sc)[1]
    if size(Sc) == (2,1)
        Sc[2,1] = α
        return
    end

    for m ∈ 2:M  # border case: first row and column
        Sc[m,1] = α * Sc[m-1,1]
        m == M ? nothing : Sc[1,m] = Sc[m,1]'
    end
    Sc[2,2] = abs2α * Sc[1,1] + Sc[1,1]
    for n ∈ 3:M # border case: second row and column
        #Sc[2,n] = Sc[1,n-1] + abs2α * Sc[1,n-1] + conj(α) * (n-2) * Sc[1,n-2]
        Sc[n, 2] = Sc[n-1,1] + abs2α * Sc[n-1, 1] + α * (n-2) * Sc[n-2, 1]
        n == M ? nothing : Sc[2,n] = Sc[n,2]'
    end
    for m ∈ 3:M
        # no conjugation for the diagonal elements
        m == M ? nothing : Sc[m,m] = Sc[m-1,m-1] + abs2α*Sc[m-1,m-1] + (m-2)*α*Sc[m-2,m-1] + (m-2)*conj(α)*Sc[m-1,m-2] + (m-2)*(m-2)*Sc[m-2,m-2]
        for n ∈ m+1:M
            #Sc[m,n] = Sc[m-1,n-1] + abs2α*Sc[m-1,n-1] + (m-2)*α*Sc[m-2,n-1] + (n-2)*conj(α)*Sc[m-1,n-2] + (m-2)*(n-2)*Sc[m-2,n-2]
            Sc[n, m] = Sc[n-1, m-1] + abs2α*Sc[n-1, m-1] + (m-2)*conj(α)*Sc[n-1, m-2] + (n-2)*α*Sc[n-2, m-1] + (m-2)*(n-2)*Sc[n-2, m-2]
            # if n == M, do not conjugate
            n == M ? nothing : Sc[m,n] = Sc[n,m]'
        end
    end
    return 
end




"""
    a!(a, S, α)

Update the matrix `a`, representing the annihilation operator with coherent field amplitude `α` and overlap matrix `S`.
"""
function a!(a::AbstractMatrix{T}, S::AbstractMatrix{T}, 
            α::T)::Nothing where {T <: Complex}

    size(a)[1] - 1 == size(a)[2] && size(a) == size(S) || error("S and a must be a matrix of size (N+1)x(N)")
    N = size(S)[2]

    a .= α * S
    a[:, 2:end] .+=  @view(S[:,1:end-1]) .* t(1:N-1)
    return
end


"""
    adag!(a, S, α)

Update the matrix `adag`, representing the creation operator with coherent field amplitude `α` and overlap matrix `S`.
"""
function adag!(adag::AbstractMatrix{T}, S::AbstractMatrix{T}, 
            α::T)::Nothing where {T <: Complex}

    size(adag)[1] - 1 == size(adag)[2] && size(adag) == size(S) || error("S and a must be a matrix of size (N+1)x(N)")
    N = size(S)[2]

    adag .= conj(α) * S
    adag[2:end, :] .+=  @view(S[1:end-1,:]) .* (1:N)
    return
end


"""
    a2!(a2, S, α)

Update the matrix `a2`, representing the operator a^2 with coherent field amplitude `α` and overlap matrix `S`.
"""
function a2!(a2::AbstractMatrix{T}, S::AbstractMatrix{T}, 
    α::T)::Nothing where {T <: Complex}

    size(a2)[1] - 1 == size(a2)[2] && size(a2) == size(S) ||error("S and a2 must be a matrix of size (N+1)x(N)")
    N = size(S)[2]

    a2 .= α^2 * S
    a2[:, 2:end] .+=  @view(S[:,1:end-1]) .* t(1:N-1) * 2*α
    if N > 2
        a2[:, 3:end] .+= @view(S[:,1:end-2]) .* t((1:N-2) .* (2:N-1))
    end
    return
end


"""
    ad_a!(ad_a, S, α)

Update the matrix `ad_a`, representing the operator a'a with coherent field amplitude `α` and overlap matrix `S`.
"""
function ad_a!(ad_a::AbstractMatrix{T}, S::AbstractMatrix{T}, 
    α::T)::Nothing where {T <: Complex}

    size(ad_a)[1] - 1 == size(ad_a)[2] && size(ad_a) == size(S)  || error("S and ad_a must be a matrix of size (N+1)x(N)")
    N = size(S)[2]

    ad_a .= abs2(α) * S
    ad_a[2:end, 2:end] .+=  @view(S[1:end-1,1:end-1]) .* (t(1:N-1) .* (1:N))

    # Maybe this can be done more efficient
    tmp1 = S .* (1:N+1) * α    
    ad_a[2:end, :] .+=  tmp1[1:end-1, :]
    tmp2 = S .* transpose(1:N) * conj(α)
    ad_a[:, 2:end] .+= tmp2[:, 1:end-1]
    return
end


"""
    ad2_a2!(ad2_a2, S, α)

Update the matrix `ad2_a2`, representing the operator a'^2 a^2 with coherent field amplitude `α` and overlap matrix `S`.
"""
function ad2_a2!(ad2_a2::AbstractMatrix{T}, S::AbstractMatrix{T}, 
    α::T)::Nothing where {T <: Complex}
    N = size(S,2)

    abs2α = abs2(α)
    ad2_a2 .= abs2α^2 * S
    ad2_a2[2:end, 2:end] .+=  @view(S[1:end-1,1:end-1]) .* (t(1:N-1) .* (1:N)) * 4*abs2α
    # tmp = @view(S[1:end-1,:]) .* (1:N) * 2*α * abs2α
    ad2_a2[2:end, :] .+=  @view(S[1:end-1,:]) .* (1:N) * 2*α * abs2α
    ad2_a2[:, 2:end] .+=  @view(S[:,1:end-1]) .* t(1:N-1) * 2*conj(α) * abs2α
    if N > 2
        ad2_a2[:,3:end] .+= @view(S[:,1:end-2]) .* t((1:N-2) .* (2:N-1)) * conj(α)^2
        ad2_a2[3:end, :] .+= @view(S[1:end-2,:]) .* (2:N) .* (1:N-1) * α^2
        ad2_a2[3:end, 2:end] .+=  @view(S[1:end-2,1:end-1]) .* t(1:N-1) .* ((1:N-1) .* (2:N)) * 2*α
        ad2_a2[2:end, 3:end] .+=  @view(S[1:end-1,1:end-2]) .*(1:N) .* t((1:N-2) .* (2:N-1)) * 2*conj(α)
        ad2_a2[3:end, 3:end] .+= @view(S[1:end-2,1:end-2]) .* t((1:N-2) .* (2:N-1)) .* (1:N-1) .* (2:N)
    end
    return
end


"""
    lmatrix!(L, B, S, H, [J, JdagJ, Jdag=adjoint.(J), rates=ones(length(J))])

Update the Liouvillian matrix with respect to the Hamiltonian `H`, density matrix `B`, overlap matrix `S` 
and jump operators `J`, `JdagJ`, `Jdag`, with the corresponding rates.

Note that L, S, H, [J], [JdagJ], [Jdag] are (n+1) x (n) matrices and B is a n x n matrix.
"""
function lmatrix!(L, B, S, H, J, JdagJ; Jdag=adjoint.(J), rates=ones(length(J)))
    size(J) == size(J) == size(Jdag) == size(rates) || error("J, Jdag, JdagJ, rates do not have the same lengths")
    size(L) == size(H) == size(S) == (size(B) .+ (1,0)) && all(==(size(H)),size.(J)) == all(==(size(H)),size.(JdagJ)) || error("dimensions of L, H, S, B, J, JdagJ do not match")

    idcs = CartesianIndices(size(B))
    Heff = copy(H) # This is important 
    @inbounds for (γ,  AdagA) ∈ zip(rates, JdagJ)
        Heff .-= 1im/2*γ*AdagA
    end

    L .= -1im * Heff * B * @view(S[idcs])
    L .+= 1im * S* B * @view(Heff[idcs])
    @inbounds for (γ, A, Adag) ∈ zip(rates, J, Jdag)
        L .+= γ * A * B * @view(Adag[idcs])
    end
    return
end


"""
    lmatrix!(L, B, S, H)

Update the Liouvillian matrix with respect to the Hamiltonian `H`, density matrix `B`, overlap matrix `S`.

Note that L, S, H are (n+1) x (n) matrices and B is a n x n matrix.
"""
function lmatrix!(L, B, S, H)
    size(L) == size(H) == size(S) == (size(B) .+ (1,0)) || error("dimensions of L, H, S, B do not match")

    idcs = CartesianIndices(size(B))
    L .= -1im * H * B * @view(S[idcs])
    L .+= 1im * S* B * @view(H[idcs])
    return
end

"""
    lmatrix_efficient!(L, B, S, H, [J, JdagJ, Jdag=adjoint.(J), rates=ones(length(J))])

(Efficiently) update the Liouvillian matrix with respect to the Hamiltonian `H`, density matrix `B`, overlap matrix `S` 
and jump operators `J`, `JdagJ`, `Jdag`, with the corresponding rates.

Note that L, S, H, [J], [JdagJ], [Jdag] are (n+1) x (n) matrices and B is a n x n matrix.
"""
function lmatrix_efficient!(L, B, S, H, J, JdagJ; Jdag=adjoint.(J), rates=ones(length(J)))
    size(J) == size(J) == size(Jdag) == size(rates) || error("J, Jdag, JdagJ, rates do not have the same lengths")
    size(L) == size(H) == size(S) == (size(B) .+ (1,0)) && all(==(size(H)),size.(J)) == all(==(size(H)),size.(JdagJ)) || error("dimensions of L, H, S, B, J, JdagJ do not match")

    idcs = CartesianIndices(size(B))
    tmp = Matrix{eltype(L)}(undef, size(B)) # temporary matrix for cache

    # Hamiltonian
    mul!(tmp, B, @view(S[idcs]), -eltype(B)(im), false) # tmp .= -1im * B * S
    mul!(L, H, tmp, true, false) # L .= H * tmp
    mul!(tmp, B, @view(H[idcs]), eltype(B)(im), false) # tmp .= 1im * B * H
    mul!(L, S, tmp, true, true) # L .= S * tmp + L

    @inbounds for i ∈ eachindex(J)
        #Quantum Jump: γᵢJᵢ*B*Jᵢ†
        mul!(tmp, B, @view(Jdag[i][idcs]), eltype(B)(rates[i]), false) # tmp .= rates[i] * B * Jdag[i]
        mul!(L, J[i], tmp, true, true) # L .= J[i] * tmp + L

        #Recycling terms: -1/2*γᵢJᵢ†*Jᵢ*B - 1/2*γᵢB*Jᵢ†*Jᵢ
        mul!(tmp, B, @view(JdagJ[i][idcs]), eltype(B)(-0.5), false) # tmp .= -1/2 * B * JdagJ[i]
        mul!(L, S, tmp, eltype(B)(rates[i]), true) # L .= γᵢ * S * tmp + L

        mul!(tmp, B, @view(S[idcs]), eltype(B)(-0.5), false) # tmp .= -1/2 * B * S
        mul!(L, JdagJ[i], tmp, eltype(B)(rates[i]), true) # L .= γᵢ * JdagJᵢ * tmp + L
    end
    return L
end


"""
    polyα(ord, α)

Return lower row of S^(left)*S⁻¹, given an order `ord` and coherent field amplitude `α`.
"""
function polyα(ord, α)
    return [binomial(ord, k-1)*α^(ord-k+1)*(-1)^k for k ∈ 1:ord]
end

function polyα!(polα, ord, α)
    @inbounds for k ∈ 1:ord
        # This could be better optimized
        polα[k] = binomial(ord, k-1)*α^(ord-k+1)*(-1)^k
    end
end


"""
    hermitian(vec[, n])

Return the Hermitian matrix of size `n` corresponding to the vector `vec`.
Here, vec contains the upper triangle with column-stacked elements.
Here n must be n = 1/2 *(sqrt(8* length(vec) + 1) - 1)
"""
function hermitian(vec::AbstractArray{T}, n::Int) where T <: Number
    #n = 1/2 *(sqrt(8* length(vec) + 1) - 1)
    k=0; 
    return Hermitian(T[ i<=j ? (k+=1; vec[k]) : 0 for i=1:n, j=1:n ])
end


function hermitian(vec::AbstractArray{T}) where T <: Number
    n = 1/2 *(sqrt(8* length(vec) + 1) - 1)
    isinteger(n) || error("Input vector does not correspond to upper triangular shape of a Hermitian matrix")
    return hermitian(vec, Int(n))
end



"""
    vec(herm)

Return the vector corresponding to the upper triangle of the Hermitian matrix `herm`.
"""
function vec(herm::Hermitian)
    return herm[triu!(trues(size(herm)))]
end