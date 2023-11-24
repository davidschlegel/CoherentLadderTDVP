mutable struct CoherentTDVP
    order::Vector{Int64}
    H::Matrix{ComplexF64} # this can be optimized
    S::Vector{Matrix{ComplexF64}}
    S⁻¹ ::Vector{Matrix{ComplexF64}}
    L::Matrix{ComplexF64}
    a:: Vector{Matrix{ComplexF64}}
    ad::Vector{Matrix{ComplexF64}}
    ad_a:: Vector{Matrix{ComplexF64}}
    a2:: Vector{Matrix{ComplexF64}}
    C0el::Vector{Int64}
    polα::Vector{Transpose{ComplexF64, Vector{ComplexF64}}}
    drift::Matrix{ComplexF64}
    Skron::Matrix{ComplexF64}    
    ad2_a2:: Vector{Matrix{ComplexF64}}
    # Add any other parameters or fields necessary

    # Constructor
    function CoherentTDVP(order::AbstractArray{T}) where T <: Int
        self = new()
        self.order = order
        self.S = [Matrix{ComplexF64}(undef, o+1, o) for o ∈ order]
        #set each first element to one.
        for i ∈ eachindex(order)
               self.S[i][1,1] = 1
        end
        self.H = Matrix{ComplexF64}(undef, prod(order .+ 1), prod(order))
        self.S⁻¹ = [Matrix{ComplexF64}(undef, o, o) for o ∈ order]
        self.L = Matrix{ComplexF64}(undef, prod(order .+ 1), prod(order))
        self.a = [Matrix{ComplexF64}(undef, o+1, o) for o ∈ order]
        self.ad = [Matrix{ComplexF64}(undef, o+1, o) for o ∈ order]        
        self.ad_a = [Matrix{ComplexF64}(undef, o+1, o) for o ∈ order]
        self.a2 = [Matrix{ComplexF64}(undef, o+1, o) for o ∈ order]
        self.ad2_a2 = [Matrix{ComplexF64}(undef, o+1, o) for o ∈ order]
        self.C0el = [factorial(o) for o ∈ order]
        self.polα = [transpose(Vector{ComplexF64}(undef, o)) for o ∈ order]
        self.drift = Matrix{ComplexF64}(undef, prod(order), prod(order))
        self.Skron = Matrix{ComplexF64}(undef, prod(order .+1), prod(order))

        # Initialize other fields as needed
        return self
    end
end



function (prob::CoherentTDVP)(du, u, p, t)

    # vectors of α
    α = u.x[1]
    # matrix B
    B = u.x[2]

    dα = du.x[1]
    dB = du.x[2]
    # update S-matrices, polyα_k, and S⁻¹
    update!(prob, α)

    # update Liouvillian and all other necessary matrices to compute L
    liouvillian!(prob, α, B, p)

    # compute derivative of alpha
    Threads.@threads for k ∈ eachindex(prob.order)
        dα[k] = alpha_derivative(prob, B, k; regnum=1e-7, regden=1e-8)
    end

    # update derivative of B
    B_derivative!(prob, dB, dα)
end


function update!(prob::CoherentTDVP, α::Vector{T}) where T
    # update S-matrixes
    for i ∈ eachindex(prob.order)
        smatrix!(prob.S[i], α[i])
        prob.S⁻¹[i] .= inv(Hermitian(@view(prob.S[i][1:end-1, :])))
        polyα!(prob.polα[i], prob.order[i], α[i])
    end
    prob.Skron .= reduce(kron, prob.S)
end    



function alpha_derivative(prob::CoherentTDVP, B::Matrix{T}, k::Int; regnum=1e-16, regden = 1e-8)::T where T

    #Numerator

    # Create a view of L of size prod(order[1:end .! = k]) x prod(order)
    # At mode k, only the last row is taken
    Lbk_row = get_rowk(prob.L, k, prob.order)

    #Polalphak_row_L = prob.polα[k]*[get_rowk(prob.L, i, prob.order) for i ∈ eachindex(prob.order)]

    # Combine into Y0k
    # Important to copy Lbk_row, otherwise it will be overwritten
    Y0k_row = copy(Lbk_row)
    # Subtract the contribution of polyalpha
    for j ∈ 1:prob.order[k]
        # Multiply scalar factor polyα_k[j] with j-th row of L associated to mode k
        Y0k_row .-=  prob.polα[k][j] * get_row_i_modek(prob.L, k, prob.order, j)
    end
    # Create a view of B of size prod(order) x prod(order[1:end .! = k])
    # At mode k, only the last column is taken
    Bk_col = get_colkb(B, k, prob.order)
    trYk = tr(Y0k_row*Bk_col)
    
    #Denominator

    #C0
    # Create a view of S of size prod(order[1:end .! = k]) x prod(order[1:end .! = k])
    # Effectively taking the Kronecker product of all S, except S[k]
    # Avoids recomputing the Kronecker product
    Skronk = get_Skronk(prob.Skron, k, prob.order)
    # Since B is Hermitian, we can conjugate transpose Bk_col' to get Bk_row
    #Bk_row = Bk_col'
    Bk_row = get_rowkb(B, k, orders) # would do the same thing
    # Make a view of S of size prod(order) x prod(order)
    S = getindexmat(prob.Skron, prob.order)
    C = Skronk*(Bk_row*S*Bk_col)
    trC = tr(C)*prob.C0el[k]
    return (trYk + regnum) / (trC + regden) 
end



function B_derivative!(prob::CoherentTDVP, dB::Matrix{T}, dα::Vector{T}) where T
    # derivative of B
    identities = [Matrix(I, o, o) for o ∈ prob.order]
    for (k,o) ∈ enumerate(prob.order)
        S_l_S⁻¹_k = diagm(o, o, (1 => ones(ComplexF64, o-1)))
        S_l_S⁻¹_k[end,:] .= transpose(prob.polα[k])
        # kron with identities, were S_l_S⁻¹_k is at position k
        # Can be optimized
        S_l_S⁻¹_k_kron = reduce(kron, [identities[1:k-1]; [S_l_S⁻¹_k]; identities[k+1:end]])
        if k == 1
            prob.drift .= conj(dα[k])*B*S_l_S⁻¹_k_kron
        else
            prob.drift .+= conj(dα[k])*B*S_l_S⁻¹_k_kron
        end
    end

    L = getindexmat(prob.L, prob.order)

    # Use `kronecker` instead?
    S⁻¹ = reduce(kron, prob.S⁻¹)

    dB .= S⁻¹*L*S⁻¹
    dB .-= prob.drift
    dB .-= prob.drift'
end