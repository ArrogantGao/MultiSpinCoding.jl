struct SpinGlass{INT}
    graph::SimpleGraph{Int}
    J::Matrix{Int}

    spins::Vector{INT}
    Js::Matrix{INT}

    n::Int

    function SpinGlass(graph::SimpleGraph{Int}, J::Matrix{Int}, n::Int)

        @assert maximum(degree(graph)) <= 6
        @assert size(J) == (nv(graph), nv(graph))
        @assert issymmetric(J)
        for e in edges(graph)
            @assert abs(J[src(e), dst(e)]) == 1
        end

        N = (n - 1) ÷ 64 + 1
        INT = BitStr{n, LongLongUInt{N}}
        spins = [random_state(n) for _ in 1:nv(graph)]
        Js = zeros(INT, nv(graph), nv(graph))
        for e in edges(graph)
            Js[src(e), dst(e)] = Js[dst(e), src(e)] = bit_J(n, J[src(e), dst(e)])
        end

        return new{INT}(graph, J, spins, Js, n)
    end
end

struct Scheduler
    nsweeps::Int
    betas::Vector{Float64}
    probs::Vector{NTuple{6, Float64}}

    function Scheduler(nsweeps::Int, beta0::Float64, beta1::Float64)
        @assert beta0 > 0
        @assert beta1 > 0
        @assert beta0 < beta1
        @assert nsweeps > 1

        betas = [beta0 + (beta1 - beta0) * i / (nsweeps - 1) for i in 0:nsweeps-1]
        probs = Vector{NTuple{6, Float64}}(undef, nsweeps)
        for (i, beta) in enumerate(betas)
            p0 = exp(-2 * beta)
            probs[i] = (p0, p0^2, p0^3, p0^4, p0^5, p0^6)
        end

        return new(nsweeps, betas, probs)
    end

    function Scheduler(betas::Vector{Float64})
        nsweeps = length(betas)
        probs = Vector{NTuple{6, Float64}}(undef, nsweeps)
        for (i, beta) in enumerate(betas)
            p0 = exp(-2 * beta)
            probs[i] = (p0, p0^2, p0^3, p0^4, p0^5, p0^6)
        end

        return new(nsweeps, betas, probs)
    end
end