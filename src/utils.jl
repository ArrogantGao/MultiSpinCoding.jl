function random_J(graph::SimpleGraph{Int})
    J = zeros(Int, nv(graph), nv(graph))
    for e in edges(graph)
        J[src(e), dst(e)] = J[dst(e), src(e)] = rand(1:2) == 1 ? 1 : -1
    end
    return J
end

# generate a LongLongUInt, n is the number of replicas
function random_state(n::Int)
    N = (n - 1) ÷ 64 + 1
    INT = BitStr{n, LongLongUInt{N}}
    spins = [rand(0:1) for _ in 1:n]
    positions = [i for i in 1:n if spins[i] == 1]
    return bmask(INT, positions...)
end

function bit_J(n::Int, j::Int)
    N = (n - 1) ÷ 64 + 1
    INT = BitStr{n, LongLongUInt{N}}
    if j == 1
        return bmask(INT, 1:n)
    elseif j == -1
        return zero(INT)
    end
end

# io
function load_instance(path::String)
    open(path, "r") do f
        n = parse(Int, readline(f))
        J = zeros(Int, n, n)
        g = SimpleGraph(n)
        for line in readlines(f)
            i, j, k = parse.(Int, split(line))
            i += 1
            j += 1
            J[i, j] = J[j, i] = k
            add_edge!(g, i, j)
        end
        return g, J
    end
end