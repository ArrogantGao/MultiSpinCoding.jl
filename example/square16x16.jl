using MultiSpinCoding

# load the graph and the coupling matrix
g, J = load_instance("square16x16.txt")

# generate the model, n = 256 is the number of replicas
model = SpinGlass(g, J, 256)

# generate the scheduler, 2000 sweeps, beta0 = 0.1, beta1 = 3.0
scheduler = Scheduler(2000, 0.1, 3.0)

# run the simulated annealing
sa!(model, scheduler)

# calculate the energies
energies = cal_energies(model)

using Plots

histogram(energies, bins=-400:1:-300, xlabel="Energy", ylabel="Frequency", title="Square 16x16")
savefig("square16x16.png")