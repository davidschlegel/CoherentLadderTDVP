function ode_function!(du::AbstractVector{T}, u::AbstractVector{T}, 
                        p::AbstractVector{N}, t::Real) where {T <: Complex, N<:Number}
    Δ = p[1]
    α = u[1]
    B = reshape(@view(u[2:end]), (order::Int, order::Int))
    # update matrices
    smatrix!(S, α)
    ad_a!(ad_a, S, α)
    
    # Your Hamiltonian here:
    H .= Δ*ad_a

    #println(α)
    # update inverse of Sc
    S⁻¹ .= inv(Hermitian(@view(S[1:end-1, :])))

    # update Liouvillian
    lmatrix!(L, B, S, H)

    # problem here: polyα seems to be wrong here
    polα = transpose(polyα(order, α))
    trY =((transpose(@view(L[end, :])) - polα*@view(L[1:end-1, :]))*@view(B[:, end]))[1]

    trC = C0_el*(transpose(@view(B[end, :]))*@view(S[1:end-1, :])*@view(B[:, end]))

    # derivative of α, including regularization
    du[1] = (trY + 1e-16) / (trC + 1e-8)

    # drift term, avoids inverse
    drift = B[:,end] * polα
    drift[:, 2:end] .+= B[:, 1:end-1] #size(B) == (1, 1) ? nothing : 
    drift .*= conj(du[1])

    # derivative of B
    du[2:end] .= reshape(S⁻¹*@view(L[1:end-1, :])*S⁻¹ - drift - drift', order^2)
end