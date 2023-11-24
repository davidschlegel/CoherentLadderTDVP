using QuantumOptics
using LinearAlgebra


"""
    BargmanBasis - Subtype of QuantumOptics.Basis
    Fields: 
        shape::Vector{Int64}
        alpha::ComplexF64

Defines a basis of Bargman states with coherent field amplitude `α` and order `ord`.

Contructor:
    `BargmanBasis(α, ord)`
"""
mutable struct BargmanBasis <: Basis
    shape::Vector{Int64}
    alpha::ComplexF64
    # Constructor
    BargmanBasis(shape::Int64, alpha) = new([shape], alpha)
    BargmanBasis(shape::Vector{Int64}, alpha) = new(shape, alpha)
end




mutable struct SingleModeCoherentTDVP
    order::Int
    H::AbstractOperator{T, T} where T <: BargmanBasis
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


mutable struct MultimodeCoherentTDVP
    order::Int64
    H::AbstractOperator
    S::Vector{AbstractOperator}
    S⁻¹::AbstractOperator
    L::AbstractOperator
    a:: AbstractOperator
    ad_a:: AbstractOperator
    a2:: AbstractOperator
    C0el::Vector{Int64}
    polα::Vector{Transpose{ComplexF64, Vector{ComplexF64}}}
    drift::AbstractOperator
end


mutable struct QOSinglemodeCoherentTDVP
    order::Int64
    H::AbstractOperator
    S::DenseOpType
    S⁻¹::DenseOpType
    L::DenseOpType
    a:: DenseOpType
    ad_a:: DenseOpType
    a2:: DenseOpType
    C0el::Int64
    polα::Transpose{ComplexF64, Vector{ComplexF64}}
    drift::DenseOpType

    function QOSinglemodeCoherentTDVP(order::Int)
        self = new()
        self.order = order
        barg_l = BargmanBasis(1.0, order+1)
        barg_r = BargmanBasis(1.0, order)
        self.S = Operator(barg_l, barg_r, Matrix{ComplexF64}(undef, order+1, order))
        self.S.data[1,1] = 1
        self.H = Operator(barg_l, barg_r, Matrix{ComplexF64}(undef, order+1, order))
        self.S⁻¹ = Operator(barg_r, barg_r,Matrix{ComplexF64}(undef, order, order))
        self.L = Operator(barg_l, barg_r,Matrix{ComplexF64}(undef, order+1, order))
        self.a = Operator(barg_l, barg_r,Matrix{ComplexF64}(undef, order+1, order))
        self.ad_a = Operator(barg_l, barg_r,Matrix{ComplexF64}(undef, order+1, order))
        self.a2 = Operator(barg_l, barg_r,Matrix{ComplexF64}(undef, order+1, order))
        self.C0el = factorial(order)
        self.polα = transpose(Vector{ComplexF64}(undef, order))
        self.drift = Operator(barg_r, barg_r,Matrix{ComplexF64}(undef, order, order))
        # Initialize other fields as needed
        return self
    end
end


function (prob::QOSinglemodeCoherentTDVP)(du, u, p, t)
    B = u
    α = B.basis_l.alpha
    # update S-matrixes
    smatrix!(prob.S.data, α)
    # update inverse
    prob.S⁻¹.data .= inv(Hermitian(@view(prob.S.data[1:end-1, :])))
    

    # update Liouvillian
    liouvillian!(prob, α, B, p)

    polyα!(prob.polα, prob.order, α)
    #prob.polα .= transpose(polyα(prob.order, α))

    trY = ((transpose(@view(prob.L.data[end, :])) - prob.polα *@view(prob.L.data[1:end-1, :]))*@view(B.data[:, end]))[1]

    trC = prob.C0el*(transpose(@view(B.data[end, :]))*@view(prob.S.data[1:end-1, :])*@view(B.data[:, end]))

    # derivative of α, including regularization
    du.x[1][1] = (trY + 1e-16) / (trC + 1e-8)

    # drift term, avoids inverse
    prob.drift = B[:,end] * prob.polα 
    prob.drift[:, 2:end] .+= B[:, 1:end-1] #size(B) == (1, 1) ? nothing : 
    prob.drift .*= conj(du.x[1][1])

    # derivative of B
    du.x[2] .= prob.S⁻¹*@view(prob.L[1:end-1, :])*prob.S⁻¹
    du.x[2] .-= prob.drift .+ prob.drift'
end




#test
α = randn(ComplexF64)
order = 10
basis_l, basis_r = BargmanBasis(order+1, α), BargmanBasis(order, α)
B_mat = zeros(ComplexF64, order, order)
B_mat[1,1] = 1
B_op = Operator(basis_r, basis_r, B_mat)


prob = QOSinglemodeCoherentTDVP(order)

B = B_op
α = B.basis_l.alpha
# update S-matrixes
smatrix!(prob.S.data, α)
# update inverse, Hermitian can be a bit faster and more stable
prob.S⁻¹.data .= inv(Hermitian(@view(prob.S.data[1:end-1, :])))


# update Liouvillian
liouvillian!(prob, α, B, p)

polyα!(prob.polα, prob.order, α)
#prob.polα .= transpose(polyα(prob.order, α))

trY = ((transpose(@view(prob.L.data[end, :])) - prob.polα *@view(prob.L.data[1:end-1, :]))*@view(B.data[:, end]))[1]

trC = prob.C0el*(transpose(@view(B.data[end, :]))*@view(prob.S.data[1:end-1, :])*@view(B.data[:, end]))

# derivative of α, including regularization
dα = (trY + 1e-16) / (trC + 1e-8)

# drift term, avoids inverse
prob.drift.data = B.data[:,end] * prob.polα 
prob.drift.data[:, 2:end] .+= B.data[:, 1:end-1] #size(B) == (1, 1) ? nothing : 
prob.drift.data .*= conj(dα)

# derivative of B
dB = prob.S⁻¹.data*@view(prob.L.data[1:end-1, :])*prob.S⁻¹.data
dB .-= (prob.drift.data .+ prob.drift.data')














α = randn(ComplexF64)
order = 10
b_left = BargmanBasis(order+1, 1.0)
b_right = BargmanBasis(order, 1.0)
B_mat = zeros(ComplexF64, order, order)
B_mat[1,1] = 1
B_op = Operator(b_right, b_right, B_mat)

b_left = BargmanBasis(10, 1.0)
b_right = BargmanBasis(9, 1.0)

mat1 = Operator(b_left, b_right, rand(ComplexF64, 10, 9))
mat2 = Operator(b_left, b_right, rand(ComplexF64, 10, 9))
comp_b_left = b_left ⊗ b_left
comp_b_right = b_right ⊗ b_right

bigmat = LazyTensor(comp_b_left, comp_b_right, [1, 2], (mat1, mat2))


#Test if DifferentialEquations works with custom density matrix type
using DifferentialEquations

function Base.similar(barg::BargmanBasis)
    return BargmanBasis(barg.shape, barg.alpha)
end

function Base.similar(op::Operator{BL, BR, T}, J=ComplexF64) where {BL<:BargmanBasis, BR<:BargmanBasis, T}
    # It is important that also the basis is reinstantiated
    return Operator(similar(op.basis_l), similar(op.basis_r), similar(op.data, J))
end


function rhs_test2_oop(du, u, p, t)
    du .= u
end


function Base.length(op::Operator{BL, BR, T}) where {BL<:BargmanBasis, BR<:BargmanBasis, T}
# It is important that also the basis is reinstantiated
    return length(op.data) + 1
end


function Base.iterate(op::Operator{BL, BR, T}, state=1) where {BL<:BargmanBasis, BR<:BargmanBasis, T}
    if state == 1
        return op.basis_l.alpha, 2
    elseif state <= length(op)
        return getindex(op.data, state-1), state+1
    else
        return nothing
    end
end


function RecursiveArrayTools.recursivefill!(op::Operator{BL, BR, T}, a) where {BL<:BargmanBasis, BR<:BargmanBasis, T}
    op.basis_l.alpha = a
    op.basis_r.alpha = a
    recursivefill!(op.data, a)
end


function Base.getindex(op::Operator{BL, BR, T}, i::CartesianIndex) where {BL<:BargmanBasis, BR<:BargmanBasis, T}
    return getindex(op.data, i)
end

function Base.getindex(op::Operator{BL, BR, T}, i::Int, j::Int) where {BL<:BargmanBasis, BR<:BargmanBasis, T}
    return getindex(op.data, i, j)
end

function Base.getindex(op::Operator{BL, BR, T}, i::Int) where {BL<:BargmanBasis, BR<:BargmanBasis, T}
    if i == 1
        return op.basis_l.alpha
    else
        return getindex(op.data, i-1)
    end
end


function Base.eachindex(op::Operator{BL, BR, T}) where {BL<:BargmanBasis, BR<:BargmanBasis, T}
    return 1:length(op)
end


function Base.keys(op::Operator{BL, BR, T}) where {BL<:BargmanBasis, BR<:BargmanBasis, T}
    #Is this dangerous?
    return eachindex(op)
end




#Broadcasting
function Base.BroadcastStyle(::Type{Operator{BL, BR, T}}) where {BL<:BargmanBasis, BR<:BargmanBasis, T}
    return Broadcast.Style{Operator{BL, BR, T}}()
end

function Base.copyto!(dest, bc::Broadcast.Broadcasted{Operator{BL, BR, T}}) where {BL<:BargmanBasis, BR<:BargmanBasis, T}
    dest.basis_l.alpha = bc.basis_l.alpha
    dest.basis_r.alpha = bc.basis_r.alpha
    copyto!(dest.data, bc.data)
end

function Base.copyto!(dest::Operator{BL, BR, T}, bc::Broadcast.Broadcasted{Nothing}) where {BL<:BargmanBasis, BR<:BargmanBasis, T}
    dest.basis_l.alpha = bc.basis_l.alpha
    dest.basis_r.alpha = bc.basis_r.alpha
    copyto!(dest.data, bc.data)
end


#tests
B_op[1,1]
#but B_op[1] gets α
B_op[1]

prob = ODEProblem(rhs_test2_oop, B_op, (0.0, 1.0))
sol   = solve(prob, RK4())

@forward typeof(B_op) (Base.length, Base.iterate, Base.getindex, Base.setindex!)


"""
Broadcasting functions from QuantumOptics
"""
# Broadcasting
Base.size(A::DataOperator) = size(A.data)
Base.size(A::DataOperator, d) = size(A.data, d)
Base.size(A::DataOperator, d::Int) = size(A.data, d) # defined to avoid method ambiguity
@inline Base.axes(A::DataOperator) = axes(A.data)
Base.broadcastable(A::DataOperator) = A

# Custom broadcasting styles
abstract type DataOperatorStyle{BL,BR} <: Broadcast.BroadcastStyle end
struct OperatorStyle{BL,BR} <: DataOperatorStyle{BL,BR} end

# Style precedence rules
Broadcast.BroadcastStyle(::Type{<:Operator{BL,BR}}) where {BL,BR} = OperatorStyle{BL,BR}()
Broadcast.BroadcastStyle(::OperatorStyle{B1,B2}, ::OperatorStyle{B3,B4}) where {B1,B2,B3,B4} = throw(IncompatibleBases())

# Out-of-place broadcasting
@inline function Base.copy(bc::Broadcast.Broadcasted{Style,Axes,F,Args}) where {BL,BR,Style<:OperatorStyle{BL,BR},Axes,F,Args<:Tuple}
    bcf = Broadcast.flatten(bc)
    bl,br = find_basis(bcf.args)
    bc_ = Broadcasted_restrict_f(bcf.f, bcf.args, axes(bcf))
    return Operator{BL,BR}(bl, br, copy(bc_))
end
find_basis(a::DataOperator, rest) = (a.basis_l, a.basis_r)

const BasicMathFunc = Union{typeof(+),typeof(-),typeof(*)}
function Broadcasted_restrict_f(f::BasicMathFunc, args::Tuple{Vararg{<:DataOperator}}, axes)
    args_ = Tuple(a.data for a=args)
    return Broadcast.Broadcasted(f, args_, axes)
end

# In-place broadcasting
@inline function Base.copyto!(dest::DataOperator{BL,BR}, bc::Broadcast.Broadcasted{Style,Axes,F,Args}) where {BL,BR,Style<:DataOperatorStyle{BL,BR},Axes,F,Args}
    axes(dest) == axes(bc) || Base.Broadcast.throwdm(axes(dest), axes(bc))
    # Performance optimization: broadcast!(identity, dest, A) is equivalent to copyto!(dest, A) if indices match
    if bc.f === identity && isa(bc.args, Tuple{<:DataOperator{BL,BR}}) # only a single input argument to broadcast!
        A = bc.args[1]
        if axes(dest) == axes(A)
            return copyto!(dest, A)
        end
    end
    # Get the underlying data fields of operators and broadcast them as arrays
    bcf = Broadcast.flatten(bc)
    bc_ = Broadcasted_restrict_f(bcf.f, bcf.args, axes(bcf))
    copyto!(dest.data, bc_)
    return dest
end
@inline Base.copyto!(A::DataOperator{BL,BR},B::DataOperator{BL,BR}) where {BL,BR} = (copyto!(A.data,B.data); A)
@inline Base.copyto!(dest::DataOperator{BL,BR}, bc::Broadcast.Broadcasted{Style,Axes,F,Args}) where {BL,BR,Style<:DataOperatorStyle,Axes,F,Args} =
    throw(IncompatibleBases())