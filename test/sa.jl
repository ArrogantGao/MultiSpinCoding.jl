using MultiSpinCoding
using Graphs, LinearAlgebra, BitBasis
using Random
using Test
Random.seed!(1234)

using MultiSpinCoding: update_d1!, update_d2!, update_d3!, update_d4!, update_d5!, update_d6!

function flip(s, nebis, Js, u, beta)
    dE = - 2 * s * dot(nebis, Js)

    if dE < 0
        return -s
    elseif u < exp(- beta * dE)
        return - s
    else
        return s
    end
end

encode(x) = x == 1 ? bit"1" : bit"0"
decode(x) = x == bit"1" ? 1 : -1

const TINT = typeof(encode(1))
beta = 0.1
p0 = exp(-2 * beta)
probs = (p0, p0^2, p0^3, p0^4, p0^5, p0^6)
us = [0.0:0.05:0.9...]

function Jmat(js)
    mat = zeros(TINT, length(js) + 1, length(js) + 1)
    for j in 1:length(js)
        mat[1, j + 1] = mat[j + 1, 1] = encode(js[j])
    end
    return mat
end

# testing all possible filps

@testset "d = 1" begin
    nebis = [2]
    for s in [-1, 1]
        for j1 in [-1, 1], s1 in [-1, 1]
            Js = Jmat([j1])
            for u in us
                sf = flip(s, [s1], [j1], u, beta)
                spins = [encode(s), encode(s1)]
                update_d1!(1, spins, Js, nebis, probs, u)
                @test spins[1] == encode(sf)
            end
        end
    end
end

@testset "d = 2" begin
    nebis = [2, 3]
    for s in [-1, 1]
        for j1 in [-1, 1], j2 in [-1, 1], s1 in [-1, 1], s2 in [-1, 1]
            Js = Jmat([j1, j2])
            for u in us
                sf = flip(s, [s1, s2], [j1, j2], u, beta)
                spins = [encode(s), encode(s1), encode(s2)]
                update_d2!(1, spins, Js, nebis, probs, u)
                @test spins[1] == encode(sf)
            end
        end
    end
end

@testset "d = 3" begin
    nebis = [2, 3, 4]
    for s in [-1, 1]
        for j1 in [-1, 1], j2 in [-1, 1], j3 in [-1, 1], s1 in [-1, 1], s2 in [-1, 1], s3 in [-1, 1]
            Js = Jmat([j1, j2, j3])
            for u in us
                sf = flip(s, [s1, s2, s3], [j1, j2, j3], u, beta)
                spins = [encode(s), encode(s1), encode(s2), encode(s3)]
                update_d3!(1, spins, Js, nebis, probs, u)
                @test spins[1] == encode(sf)
            end
        end
    end
end

@testset "d = 4" begin
    nebis = [2, 3, 4, 5]
    for s in [-1, 1]
        for j1 in [-1, 1], j2 in [-1, 1], j3 in [-1, 1], j4 in [-1, 1], s1 in [-1, 1], s2 in [-1, 1], s3 in [-1, 1], s4 in [-1, 1]
            Js = Jmat([j1, j2, j3, j4])
            for u in us
                sf = flip(s, [s1, s2, s3, s4], [j1, j2, j3, j4], u, beta)
                spins = [encode(s), encode(s1), encode(s2), encode(s3), encode(s4)]
                update_d4!(1, spins, Js, nebis, probs, u)
                @test spins[1] == encode(sf)
            end
        end
    end
end

@testset "d = 5" begin
    nebis = [2, 3, 4, 5, 6]
    for s in [-1, 1]
        for j1 in [-1, 1], j2 in [-1, 1], j3 in [-1, 1], j4 in [-1, 1], j5 in [-1, 1], s1 in [-1, 1], s2 in [-1, 1], s3 in [-1, 1], s4 in [-1, 1], s5 in [-1, 1]
            Js = Jmat([j1, j2, j3, j4, j5])
            for u in us
                sf = flip(s, [s1, s2, s3, s4, s5], [j1, j2, j3, j4, j5], u, beta)
                spins = [encode(s), encode(s1), encode(s2), encode(s3), encode(s4), encode(s5)]
                update_d5!(1, spins, Js, nebis, probs, u)
                @test spins[1] == encode(sf)
            end
        end
    end
end

@testset "d = 6" begin
    nebis = [2, 3, 4, 5, 6, 7]
    for s in [-1, 1]
        for j1 in [-1, 1], j2 in [-1, 1], j3 in [-1, 1], j4 in [-1, 1], j5 in [-1, 1], j6 in [-1, 1], s1 in [-1, 1], s2 in [-1, 1], s3 in [-1, 1], s4 in [-1, 1], s5 in [-1, 1], s6 in [-1, 1]
            Js = Jmat([j1, j2, j3, j4, j5, j6])
            for u in us
                sf = flip(s, [s1, s2, s3, s4, s5, s6], [j1, j2, j3, j4, j5, j6], u, beta)
                spins = [encode(s), encode(s1), encode(s2), encode(s3), encode(s4), encode(s5), encode(s6)]
                update_d6!(1, spins, Js, nebis, probs, u)
                @test spins[1] == encode(sf)
            end
        end
    end
end

@testset "sa" begin
    g, J = load_instance("../example/square16x16.txt")
    model = SpinGlass(g, J, 256)
    scheduler = Scheduler(2000, 0.1, 3.0)
    sa!(model, scheduler)
    energies = cal_energies(model)
    @test minimum(energies) == -360
end