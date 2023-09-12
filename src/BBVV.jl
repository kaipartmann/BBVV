module BBVV

using Printf
using ProgressMeter
using WriteVTK
using Base.Threads
using TimerOutputs
using ThreadPinning

export PointCloud, BondBasedMaterial, VelocityBC, PreCrack, simulation

pinthreads(:numa; force=false)

const TO = TimerOutput()

struct PointCloud
    n_points::Int
    position::Matrix{Float64}
    volume::Vector{Float64}
end

struct PreCrack
    point_id_set_a::Vector{Int}
    point_id_set_b::Vector{Int}
end

struct BondBasedMaterial
    δ::Float64
    bc::Float64
    rho::Float64
    εc::Float64
end

struct VelocityBC
    fun::Function
    point_id_set::Vector{Int}
    dim::Int
end

struct SimulationParameters
    n_threads::Int
    mat::BondBasedMaterial
    pc::PointCloud
    cells::Vector{MeshCell{VTKCellType, Tuple{Int64}}}
    n_bonds::Int
    bonds::Vector{Int}
    init_dists::Vector{Float64}
    n_family_members::Vector{Int}
    hood_range::Vector{UnitRange{Int}}
    bcs::Vector{VelocityBC}
    n_timesteps::Int
    export_freq::Int
    export_path::String
end

struct GlobalStorage
    position::Matrix{Float64}
    displacement::Matrix{Float64}
    velocity::Matrix{Float64}
    velocity_half::Matrix{Float64}
    acceleration::Matrix{Float64}
    bond_failure::Vector{Int}
    damage::Vector{Float64}
    b_int::Matrix{Float64}
    n_active_family_members::Vector{Int}
end

struct ThreadLocalStorage
    tl_points::UnitRange{Int}
    pointmap::Dict{Int,Int}
    b_int::Matrix{Float64}
    n_active_family_members::Vector{Int}
end


"""
MAIN FUNCTION
"""
function simulation(pc::PointCloud, mat::BondBasedMaterial, bcs::Vector{VelocityBC};
                    precracks::Vector{PreCrack}=Vector{PreCrack}(), n_timesteps::Int=1000,
                    export_freq::Int=10, export_path::String="results")
    println("----------------------------")
    println("  PERIDYNAMICS SIMULATION")
    println("----------------------------")
    reset_timer!(TO)
    walltime = @elapsed begin
        print("initialization...")
        sp, gs, vtls = init_simulation(pc, mat, bcs, precracks, n_timesteps, export_freq,
                                       export_path)
        println("\r✔ initialization   ")
        print("time loop...")
        time_loop!(gs, vtls, sp)
        println("\r✔ time loop   ")
    end
    open(joinpath(export_path,"logfile.log"), "w+") do io
        write(io, "simulation completed after $walltime seconds (wall time)\n\n")
        show(IOContext(io, :displaysize => (24,150)), TO)
    end
    println("--- SIMULATION COMPLETED ---")
end


@timeit TO function init_simulation(pc::PointCloud, mat::BondBasedMaterial,
                                    bcs::Vector{VelocityBC}, precracks::Vector{PreCrack},
                                    n_timesteps::Int, export_freq::Int, export_path::String)
    n_threads = nthreads()

    # create SimulationParameters
    point_dist = defaultdist(pc.n_points, n_threads)
    bonds, n_family_members, init_dists = find_bonds(pc, mat.δ, point_dist, n_threads)
    n_bonds = length(bonds)
    cells = get_cells(pc.n_points)
    hood_range = find_hood_range(n_family_members, pc.n_points)
    sp = SimulationParameters(n_threads, mat, pc, cells, n_bonds, bonds, init_dists,
                              n_family_members, hood_range, bcs, n_timesteps, export_freq,
                              export_path)

    # create Vector{ThreadLocalStorage}
    vtls = Vector{ThreadLocalStorage}(undef, n_threads)
    @threads for tid in 1:n_threads
        tl_points = point_dist[tid]
        pointmap = Dict{Int,Int}()
        n_local_points = 0
        for i in tl_points
            if i in keys(pointmap)
                error("key $i is already defined in the pointmap dict!\n")
            end
            n_local_points += 1
            pointmap[i] = n_local_points
        end
        @assert length(pointmap) == n_local_points
        b_int_local = zeros(3, n_local_points)
        n_active_family_members_local = zeros(n_local_points)
        for (gi, li) in pointmap
            if li > n_local_points
                error("local_index > n_local_points!\n")
            end
            n_active_family_members_local[li] = n_family_members[gi]
        end
        vtls[tid] = ThreadLocalStorage(tl_points, pointmap, b_int_local,
                                       n_active_family_members_local)
    end

    # create GlobalStorage
    position = copy(pc.position)
    displacement = zeros(3, pc.n_points)
    velocity = zeros(3, pc.n_points)
    velocity_half = zeros(3, pc.n_points)
    acceleration = zeros(3, pc.n_points)
    b_int = zeros(3, pc.n_points)
    bond_failure = ones(Int, n_bonds)
    damage = zeros(pc.n_points)
    n_active_family_members = copy(n_family_members)
    gs = GlobalStorage(position, displacement, velocity, velocity_half, acceleration,
                       bond_failure, damage, b_int,n_active_family_members)

    define_precracks!(gs, vtls, sp, precracks)

    return sp, gs, vtls
end

@timeit TO function time_loop!(gs::GlobalStorage, vtls::Vector{ThreadLocalStorage},
                               sp::SimulationParameters)
    Δt = calc_stable_timestep(sp, vtls)
    Δt½ = 0.5Δt
    export_vtk(sp, gs, 0, 0.0)

    for t in 1:sp.n_timesteps
        time = t * Δt
        update_velocity_half!(gs, sp, Δt½)
        apply_bcs!(gs, sp, time)
        update_disp_and_position!(gs, sp, Δt)
        compute_forcedensity!(vtls, gs, sp)
        reduce_tls_to_gs!(gs, vtls)
        calc_damage!(gs, sp)
        compute_equation_of_motion!(gs, sp, Δt½)
        if mod(t, sp.export_freq) == 0
            export_vtk(sp, gs, t, time)
        end
    end

    return nothing
end

function PointCloud(lx::Real, ly::Real, lz::Real, Δx::Real)
    _gridx = range(; start = (-lx + Δx) / 2, stop = (lx - Δx) / 2, step = Δx)
    gridx = _gridx .- sum(_gridx) / length(_gridx)
    _gridy = range(; start = (-ly + Δx) / 2, stop = (ly - Δx) / 2, step = Δx)
    gridy = _gridy .- sum(_gridy) / length(_gridy)
    _gridz = range(; start = (-lz + Δx) / 2, stop = (lz - Δx) / 2, step = Δx)
    gridz = _gridz .- sum(_gridz) / length(_gridz)
    _position = vec(collect(Iterators.product(gridx, gridy, gridz)))
    position = reinterpret(reshape, eltype(eltype(_position)), _position)
    if isempty(position)
        err_msg = "Size of Δx too big to create point cloud! Decrease the value!\n"
        throw(ArgumentError(err_msg))
    end
    n_points = size(position, 2)
    volume = fill(Δx^3, n_points)
    return PointCloud(n_points, position, volume)
end

function defaultdist(sz::Int, nc::Int)
    if sz >= nc
        chunk_size = div(sz, nc)
        remainder = rem(sz, nc)
        sidx = zeros(Int64, nc + 1)
        for i in 1:(nc + 1)
            sidx[i] += (i - 1) * chunk_size + 1
            if i <= remainder
                sidx[i] += i - 1
            else
                sidx[i] += remainder
            end
        end
        grid = fill(0:0, nc)
        for i in 1:nc
            grid[i] = sidx[i]:(sidx[i + 1] - 1)
        end
        return grid
    else
        sidx = [1:(sz + 1);]
        grid = fill(0:-1, nc)
        for i in 1:sz
            grid[i] = sidx[i]:(sidx[i + 1] - 1)
        end
        return grid
    end
end

@timeit TO function find_bonds(pc::PointCloud, δ::Float64,
                               owned_points::Vector{UnitRange{Int}}, n_threads::Int)
    _bonds = Vector{Vector{Int}}(undef, n_threads)
    _init_dists = Vector{Vector{Float64}}(undef, n_threads)
    n_family_members = zeros(Int, pc.n_points)
    @threads :static for tid in 1:n_threads
        local_bonds = Vector{Int}()
        local_init_dists = Vector{Float64}()
        for i in owned_points[tid]
            num = 0
            for j in 1:pc.n_points
                if i !== j
                    idst = sqrt((pc.position[1, j] - pc.position[1, i])^2 +
                                (pc.position[2, j] - pc.position[2, i])^2 +
                                (pc.position[3, j] - pc.position[3, i])^2)
                    if idst <= δ
                        num += 1
                        push!(local_bonds, j)
                        push!(local_init_dists, idst)
                    end
                end
            end
            n_family_members[i] = num
        end
        _bonds[tid] = local_bonds
        _init_dists[tid] = local_init_dists
    end
    bonds = reduce(append!, _bonds)
    init_dists = reduce(append!, _init_dists)
    return bonds, n_family_members, init_dists
end

get_cells(n::Int) = [MeshCell(VTKCellTypes.VTK_VERTEX, (i,)) for i in 1:n]

@timeit TO function calc_stable_timestep(sp::SimulationParameters, tls::Vector{ThreadLocalStorage})
    _Δt = zeros(Float64, sp.n_threads)
    @threads :static for tid in 1:sp.n_threads
        timesteps = fill(typemax(Float64), sp.pc.n_points)
        for i in tls[tid].tl_points
            dtsum = 0.0
            for nid in sp.hood_range[i]
                j = sp.bonds[nid]
                L = sp.init_dists[nid]
                dtsum += sp.pc.volume[j] * sp.mat.bc / L
            end
            timesteps[i] = sqrt(2 * sp.mat.rho / dtsum)
        end
        _Δt[tid] = 0.7 * minimum(timesteps)
    end
    Δt = minimum(_Δt)
    return Δt
end

@timeit TO function export_vtk(sp::SimulationParameters, gs::GlobalStorage, timestep::Int,
                               time::Float64)
    filename = joinpath(sp.export_path, @sprintf("timestep_%04d", timestep))
    vtk_grid(filename, sp.pc.position, sp.cells) do vtk
        vtk["Damage", VTKPointData()] = gs.damage
        vtk["Displacement", VTKPointData()] = gs.displacement
        vtk["Time", VTKFieldData()] = time
    end
    return nothing
end

function find_hood_range(n_family_members::Vector{Int}, n_points::Int)
    hood_range = fill(0:0, n_points)
    cbond = 1
    for i in 1:n_points
        hood_range[i] = cbond:(cbond + n_family_members[i] - 1)
        cbond += n_family_members[i]
    end
    return hood_range
end

@timeit TO function define_precracks!(gs::GlobalStorage, vtls::Vector{ThreadLocalStorage},
                                      sp::SimulationParameters, precracks::Vector{PreCrack})
    for precrack in precracks
        @threads :static for tls in vtls
            tls.n_active_family_members .= 0
            for i in tls.tl_points
                for cbond in sp.hood_range[i]
                    j = sp.bonds[cbond]
                    i_is_in_set_a = in(i, precrack.point_id_set_a)
                    i_is_in_set_b = in(i, precrack.point_id_set_b)
                    j_is_in_set_a = in(j, precrack.point_id_set_a)
                    j_is_in_set_b = in(j, precrack.point_id_set_b)
                    if i_is_in_set_a && j_is_in_set_b || i_is_in_set_b && j_is_in_set_a
                        gs.bond_failure[cbond] = 0
                    end
                    li = tls.pointmap[i]
                    tls.n_active_family_members[li] += gs.bond_failure[cbond]
                end
            end
        end
    end

    # reduction necessary
    gs.n_active_family_members .= 0
    for tls in vtls
        for (gi, li) in tls.pointmap
            gs.n_active_family_members[gi] = tls.n_active_family_members[li]
        end
    end

    # initial calc_damage call
    calc_damage!(gs, sp)

    return nothing
end

@timeit TO function compute_forcedensity!(vtls::Vector{ThreadLocalStorage},
                                          gs::GlobalStorage, sp::SimulationParameters)
    @threads :static for tls in vtls
        tls.b_int .= 0
        tls.n_active_family_members .= 0
        for i in tls.tl_points
            for cbond in sp.hood_range[i]
                j = sp.bonds[cbond]
                L = sp.init_dists[cbond]
                li = tls.pointmap[i]
                ϑx = gs.position[1, j] - gs.position[1, i]
                ϑy = gs.position[2, j] - gs.position[2, i]
                ϑz = gs.position[3, j] - gs.position[3, i]
                l = sqrt(ϑx * ϑx + ϑy * ϑy + ϑz * ϑz)
                ε = (l - L) / L
                temp = gs.bond_failure[cbond] * sp.mat.bc * ε / l * sp.pc.volume[j]
                tls.b_int[1, li] += temp * ϑx
                tls.b_int[2, li] += temp * ϑy
                tls.b_int[3, li] += temp * ϑz

                # failure mechanism
                if ε > sp.mat.εc
                    gs.bond_failure[cbond] = 0
                end
                tls.n_active_family_members[li] += gs.bond_failure[cbond]
            end
        end
    end
    return nothing
end

@timeit TO function apply_bcs!(gs::GlobalStorage, sp::SimulationParameters, time::Float64)
    @sync for bc in sp.bcs
        Threads.@spawn begin
            value = bc.fun(time)
            dim = bc.dim
            @inbounds @simd for i in bc.point_id_set
                gs.velocity_half[dim, i] = value
            end
        end
    end
    return nothing
end

@timeit TO function reduce_tls_to_gs!(gs::GlobalStorage, vtls::Vector{ThreadLocalStorage})
    # reduction
    gs.b_int .= 0
    gs.n_active_family_members .= 0
    for tls in vtls
        for (gi, li) in tls.pointmap
            for d in 1:3
                gs.b_int[d, gi] = tls.b_int[d, li]
            end
            gs.n_active_family_members[gi] = tls.n_active_family_members[li]
        end
    end
    return nothing
end

@timeit TO function update_disp_and_position!(gs::GlobalStorage, sp::SimulationParameters,
                                              Δt::Float64)
    @threads :static for i in 1:sp.pc.n_points
        for d in 1:3
            gs.displacement[d, i] += gs.velocity_half[d, i] * Δt
            gs.position[d, i] += gs.velocity_half[d, i] * Δt
        end
    end
    return nothing
end

@timeit TO function update_velocity_half!(gs::GlobalStorage, sp::SimulationParameters,
                                          Δt½::Float64)
    @threads :static for i in 1:sp.pc.n_points
        for d in 1:3
            gs.velocity_half[d, i] = gs.velocity[d, i] + gs.acceleration[d, i] * Δt½
        end
    end
    return nothing
end

@timeit TO function compute_equation_of_motion!(gs::GlobalStorage, sp::SimulationParameters,
                                                Δt½::Float64)
    @threads :static for i in 1:sp.pc.n_points
        for d in 1:3
            gs.acceleration[d, i] = (gs.b_int[d, i]) / sp.mat.rho
            gs.velocity[d, i] = gs.velocity_half[d, i] + gs.acceleration[d, i] * Δt½
        end
    end
    return nothing
end

@timeit TO function calc_damage!(gs::GlobalStorage, sp::SimulationParameters)
    @threads :static for i in 1:sp.pc.n_points
        gs.damage[i] = 1 - gs.n_active_family_members[i] / sp.n_family_members[i]
    end
    return nothing
end

end
