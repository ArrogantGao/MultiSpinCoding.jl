@inbounds function sa!(model::SpinGlass, scheduler::Scheduler)
    for nsweep in 1:scheduler.nsweeps
        # nsweep % 100 == 0 && println("nsweep: $nsweep")
        for i in 1:nv(model.graph)
            # update the j-th spin
            d = degree(model.graph, i)
            spins = model.spins
            Js = model.Js
            nebis = neighbors(model.graph, i)
            probs = scheduler.probs[nsweep]
            u = rand()
            if d == 1
                update_d1!(i, spins, Js, nebis, probs, u)
            elseif d == 2
                update_d2!(i, spins, Js, nebis, probs, u)
            elseif d == 3
                update_d3!(i, spins, Js, nebis, probs, u)
            elseif d == 4
                update_d4!(i, spins, Js, nebis, probs, u)
            elseif d == 5
                update_d5!(i, spins, Js, nebis, probs, u)
            elseif d == 6
                update_d6!(i, spins, Js, nebis, probs, u)
            else
                error("degree of site $i is $d larger than 6, which is not supported")
            end
        end
    end
end

@inbounds function update_d1!(i::Int, spins::Vector{INT}, Js::Matrix{INT}, nebis::Vector{Int}, probs::NTuple{6, Float64}, u::Float64) where INT
    j = nebis[1]
    l1 = Js[i, j] ⊻ (spins[i] ⊻ spins[j])

    if u > probs[1]
        # filp the site with energy decrease
        spins[i] = spins[i] ⊻ l1
    else
        # filp all
        spins[i] = neg(spins[i])
    end
    nothing
end

@inbounds function update_d2!(i::Int, spins::Vector{INT}, Js::Matrix{INT}, nebis::Vector{Int}, probs::NTuple{6, Float64}, u::Float64) where INT
    l0 = Js[i, nebis[1]] ⊻ (spins[i] ⊻ spins[nebis[1]])
    l1 = Js[i, nebis[2]] ⊻ (spins[i] ⊻ spins[nebis[2]])

    j0 = l0 ⊻ l1
    j1 = l0 & l1

    if u > probs[2]
        mask = j0 | j1
        spins[i] = spins[i] ⊻ mask
    else
        spins[i] = neg(spins[i])
    end
    nothing
end

@inbounds function update_d3!(i::Int, spins::Vector{INT}, Js::Matrix{INT}, nebis::Vector{Int}, probs::NTuple{6, Float64}, u::Float64) where INT
    l0 = Js[i, nebis[1]] ⊻ (spins[i] ⊻ spins[nebis[1]])
    l1 = Js[i, nebis[2]] ⊻ (spins[i] ⊻ spins[nebis[2]])
    l2 = Js[i, nebis[3]] ⊻ (spins[i] ⊻ spins[nebis[3]])

    j1 = l0 ⊻ l1
    j0 = j1 ⊻ l2
    j1 = (l0 & l1) ⊻ (j1 & l2)

    if u > probs[1]
        mask = j1
        spins[i] = spins[i] ⊻ mask
    elseif u > probs[3]
        mask = j1 | j0
        spins[i] = spins[i] ⊻ mask
    else
        spins[i] = neg(spins[i])
    end
    nothing
end

@inbounds function update_d4!(i::Int, spins::Vector{INT}, Js::Matrix{INT}, nebis::Vector{Int}, probs::NTuple{6, Float64}, u::Float64) where INT
    l0 = Js[i, nebis[1]] ⊻ (spins[i] ⊻ spins[nebis[1]])
    l1 = Js[i, nebis[2]] ⊻ (spins[i] ⊻ spins[nebis[2]])
    l2 = Js[i, nebis[3]] ⊻ (spins[i] ⊻ spins[nebis[3]])
    l3 = Js[i, nebis[4]] ⊻ (spins[i] ⊻ spins[nebis[4]])

    j0 = l0 ⊻ l1
    j1 = l0 & l1
    j2 = l2 ⊻ l3
    j3 = l2 & l3

    if u > probs[2]
        mask = j1 | j3 | (j0 & j2)
        spins[i] = spins[i] ⊻ mask
    elseif u > probs[4]
        mask = j1 | j3 | j0 | j2
        spins[i] = spins[i] ⊻ mask
    else
        spins[i] = neg(spins[i])
    end
    nothing
end

@inbounds function update_d5!(i::Int, spins::Vector{INT}, Js::Matrix{INT}, nebis::Vector{Int}, probs::NTuple{6, Float64}, u::Float64) where INT
    l0 = Js[i, nebis[1]] ⊻ (spins[i] ⊻ spins[nebis[1]])
    l1 = Js[i, nebis[2]] ⊻ (spins[i] ⊻ spins[nebis[2]])
    l2 = Js[i, nebis[3]] ⊻ (spins[i] ⊻ spins[nebis[3]])
    l3 = Js[i, nebis[4]] ⊻ (spins[i] ⊻ spins[nebis[4]])
    l4 = Js[i, nebis[5]] ⊻ (spins[i] ⊻ spins[nebis[5]])

    j1 = l0 ⊻ l1
    j0 = j1 ⊻ l2
    j1 = (l0 & l1) ⊻ (j1 & l2)
    j2 = l3 ⊻ l4
    j3 = l3 & l4

    if u > probs[1]
        mask = ((j1 | j3) & (j0 | j2)) | (j1 & j3)
        spins[i] = spins[i] ⊻ mask
    elseif u > probs[3]
        mask = (j0 & j2) | j1 | j3
        spins[i] = spins[i] ⊻ mask
    elseif u > probs[5]
        mask = j0 | j2 | j1 | j3
        spins[i] = spins[i] ⊻ mask
    else
        spins[i] = neg(spins[i])
    end
    nothing    
end

@inbounds function update_d6!(i::Int, spins::Vector{INT}, Js::Matrix{INT}, nebis::Vector{Int}, probs::NTuple{6, Float64}, u::Float64) where INT
    l0 = Js[i, nebis[1]] ⊻ (spins[i] ⊻ spins[nebis[1]])
    l1 = Js[i, nebis[2]] ⊻ (spins[i] ⊻ spins[nebis[2]])
    l2 = Js[i, nebis[3]] ⊻ (spins[i] ⊻ spins[nebis[3]])
    l3 = Js[i, nebis[4]] ⊻ (spins[i] ⊻ spins[nebis[4]])
    l4 = Js[i, nebis[5]] ⊻ (spins[i] ⊻ spins[nebis[5]])
    l5 = Js[i, nebis[6]] ⊻ (spins[i] ⊻ spins[nebis[6]])

    j1 = l0 ⊻ l1
    j0 = j1 ⊻ l2
    j1 = (l0 & l1) ⊻ (j1 & l2)

    j3 = l3 ⊻ l4
    j2 = j3 ⊻ l5
    j3 = (l3 & l4) ⊻ (j3 & l5)

    if u > probs[2]
        mask = ((j1 | j3) & (j0 | j2)) | (j1 & j3)
        spins[i] = spins[i] ⊻ mask
    elseif u > probs[4]
        mask = (j0 & j2) | j1 | j3
        spins[i] = spins[i] ⊻ mask
    elseif u > probs[6]
        mask = j0 | j2 | j1 | j3
        spins[i] = spins[i] ⊻ mask
    else
        spins[i] = neg(spins[i])
    end
    nothing
end

function cal_energies(model::SpinGlass)
    spins = model.spins
    n = model.n
    g = model.graph
    J = model.J

    energies = zeros(Float64, n)
    for e in edges(g)
        for i in 1:n
            j = J[src(e), dst(e)]
            spin1 = spins[src(e)]
            spin2 = spins[dst(e)]
            si = readbit(spin1, i)
            sj = readbit(spin2, i)
            sisj = si == sj ? 1 : -1
            energies[i] += j * sisj
        end
    end
    return energies
end