using BBVV

RESPATH::String = length(ARGS) ≥ 1 ? ARGS[1] : joinpath(@__DIR__, "results")

function bbvv_benchmark(simname, l, Δx, v0, nt)
    pc = PointCloud(l, l, 0.1l, Δx)
    δ = 3.015Δx
    E = 2.1e5
    nu = 0.25
    bc = 18 * E / (3 * (1 - 2 * nu)) / (π * δ^4)
    rho = 8e-6
    εc = 0.01
    mat = BondBasedMaterial(δ, bc, rho, εc)
    a = 0.5l
    set_a = findall(p -> p[1] ≤ -l/2+a && 0 ≤ p[2] ≤ 2δ, eachcol(pc.position))
    set_b = findall(p -> p[1] ≤ -l/2+a && -2δ ≤ p[2] < 0, eachcol(pc.position))
    precracks = [PreCrack(set_a, set_b)]
    set_top = findall(p -> p[2] ≥ l/2-Δx, eachcol(pc.position))
    set_bottom = findall(p -> p[2] ≤ -l/2+Δx, eachcol(pc.position))
    bc_bottom = VelocityBC(t -> -v0, set_bottom, 2)
    bc_top = VelocityBC(t -> v0, set_top, 2)
    bcs = [bc_bottom, bc_top]
    path = joinpath(RESPATH, simname)
    ispath(path) && rm(path; recursive=true, force=true)
    !ispath(path) && mkpath(path) # create the path if it does not exist
    simulation(pc, mat, bcs; precracks=precracks, n_timesteps=nt, export_freq=10000,
               export_path=path)
    return nothing
end

##--
simname = "bbvv_" * string(Threads.nthreads())
l = 1.0
Δx = 1/20
v0 = 10
nt = 1000
println("\n** compilationrun $simname")
@time bbvv_benchmark(simname, l, Δx, v0, nt) # compilation run

##--
Δx = 1/150
nt = 2000
println("\n** benchmarkrun $simname")
@time bbvv_benchmark(simname, l, Δx, v0, nt)
