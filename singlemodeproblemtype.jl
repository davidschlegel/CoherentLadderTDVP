# Here, we define the problems including initialization processes

mutable struct SingleModeCoherentTDVP
    order::Int
    H::Matrix{ComplexF64}
    S::Matrix{ComplexF64}
    S⁻¹ ::Matrix{ComplexF64}
    L::Matrix{ComplexF64}
    a:: Matrix{ComplexF64}
    ad_a:: Matrix{ComplexF64}
    a2:: Matrix{ComplexF64}
    C0el::Int
    polα::Transpose{ComplexF64, Vector{ComplexF64}}
    drift::Matrix{ComplexF64}
    #ad2_a2:: Matrix{ComplexF64}
    # Add any other parameters or fields necessary

    # Constructor
    function SingleModeCoherentTDVP(order::Int)
        self = new()
        self.order = order
        self.S = Matrix{ComplexF64}(undef, order+1, order)
        self.S[1,1] = 1
        self.H = Matrix{ComplexF64}(undef, order+1, order)
        self.S⁻¹ = Matrix{ComplexF64}(undef, order, order)
        self.L = Matrix{ComplexF64}(undef, order+1, order)
        self.a = Matrix{ComplexF64}(undef, order+1, order)
        self.ad_a = Matrix{ComplexF64}(undef, order+1, order)
        self.a2 = Matrix{ComplexF64}(undef, order+1, order)
        self.C0el = factorial(order)
        self.polα = transpose(Vector{ComplexF64}(undef, order))
        self.drift = Matrix{ComplexF64}(undef, order, order)
        # Initialize other fields as needed
        return self
    end
end



function (prob::SingleModeCoherentTDVP)(du, u, p, t)
    α = u.x[1][1]
    B = u.x[2]
    # update S-matrixes
    smatrix!(prob.S, α)
    prob.S⁻¹ .= inv(Hermitian(@view(prob.S[1:end-1, :])))
    

    # update Liouvillian
    liouvillian!(prob, α, B, p)

    polyα!(prob.polα, prob.order, α)
    #prob.polα .= transpose(polyα(prob.order, α))

    trY = ((transpose(@view(prob.L[end, :])) - prob.polα *@view(prob.L[1:end-1, :]))*@view(B[:, end]))[1]

    trC = prob.C0el*(transpose(@view(B[end, :]))*@view(prob.S[1:end-1, :])*@view(B[:, end]))

    # derivative of α, including regularization
    du.x[1][1] = (trY + 1e-16) / (trC + 1e-8)

    # drift term, avoids inverse
    prob.drift = B[:,end] * prob.polα 
    prob.drift[:, 2:end] .+= B[:, 1:end-1] #size(B) == (1, 1) ? nothing : 
    prob.drift .*= conj(du.x[1][1])

    # derivative of B
    du.x[2] .= prob.S⁻¹*@view(prob.L[1:end-1, :])*prob.S⁻¹
    du.x[2] .-= prob.drift
    du.x[2] .-= prob.drift'
end