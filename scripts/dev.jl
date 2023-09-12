using BBVV
using Base.Threads
using Printf

##--
RESPATH = joinpath(@__DIR__, "results")
simname = "run1"
l = 1.0
Δx = 1/60
v0 = 10
n_timesteps = 1000
export_freq = 10
export_path = joinpath(RESPATH, simname)

##--
pc = PointCloud(l, l, 0.1l, Δx)
δ = 3.015Δx
E = 2.1e5
nu = 0.25
bc = 18 * E / (3 * (1 - 2 * nu)) / (π * δ^4)
rho = 8e-6
εc = 0.005
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
!ispath(path) && mkpath(path) # create the path if it doe

##--
sp, gs, vtls = BBVV.init_simulation(pc, mat, bcs, precracks, n_timesteps, export_freq, export_path)

##--
function workload_dist(sp, vtls)
    n_threads = length(vtls)
    n_bonds_per_thread = zeros(Int, n_threads)
    @threads :static for tid in 1:n_threads
        n_bonds_per_thread[tid] = sum(sp.n_family_members[vtls[tid].tl_points])
    end
    return n_bonds_per_thread
end

@time wd = workload_dist(sp, vtls)
wd ./ maximum(wd)


##--
function workload_of_dist(point_distribution::Vector{UnitRange{Int}},
                          n_family_members::Vector{Int})
    n_threads = length(point_distribution)
    n_bonds_per_thread = zeros(Int, n_threads)
    @threads :static for tid in 1:n_threads
        thread_local_points = point_distribution[tid]
        n_bonds_per_thread[tid] = sum(n_family_members[thread_local_points])
    end
    wl_pct = n_bonds_per_thread ./ maximum(n_bonds_per_thread)
    println("thread | number of bonds | percentage")
    for tid in 1:n_threads
        @printf("%6d | %15d | %8.2f %%\n", tid, n_bonds_per_thread[tid], 100wl_pct[tid])
    end
    return n_bonds_per_thread
end
workload_of_dist(point_ddist, n_family_members);

##-- balancedist - script
# n_family_members = copy(sp.n_family_members)
# n_threads = nthreads()
# n_points = length(n_family_members)
# n_bonds = sum(n_family_members)
# n_bonds_per_thread_target = n_bonds ÷ n_threads
# point_ddist = BBVV.defaultdist(n_points, n_threads)
# # n_bonds_per_thread = workload_of_dist(point_ddist, n_family_members)
# new_dist = fill(0:0, n_threads)
# start_point = 1
# for tid in 1:n_threads
#     end_point = last(point_ddist[tid])
#     n_bonds_per_tid = sum(n_family_members[start_point:end_point])
#     if n_bonds_per_tid < n_bonds_per_thread_target
#         while n_bonds_per_tid < n_bonds_per_thread_target && end_point < n_points
#             end_point += 1
#             n_bonds_per_tid = sum(n_family_members[start_point:end_point])
#         end
#     elseif n_bonds_per_tid > n_bonds_per_thread_target
#         while n_bonds_per_tid > n_bonds_per_thread_target && end_point < n_points
#             end_point -= 1
#             n_bonds_per_tid = sum(n_family_members[start_point:end_point])
#         end
#     end
#     new_dist[tid] = start_point:end_point
#     start_point = end_point + 1
# end
# workload_of_dist(new_dist, n_family_members)

##--
function balancedist(point_ddist, n_family_members)
    n_threads = length(point_ddist)
    n_points = length(n_family_members)
    n_bonds = sum(n_family_members)
    n_bonds_per_thread_target = n_bonds ÷ n_threads
    balancedist = fill(0:0, n_threads)
    start_point = 1
    for tid in 1:n_threads
        end_point = last(point_ddist[tid])
        n_bonds_per_tid = sum(n_family_members[start_point:end_point])
        if n_bonds_per_tid < n_bonds_per_thread_target && end_point < n_points
            while n_bonds_per_tid < n_bonds_per_thread_target
                end_point += 1
                n_bonds_per_tid += n_family_members[end_point]
            end
        elseif n_bonds_per_tid > n_bonds_per_thread_target && end_point < n_points
            while n_bonds_per_tid > n_bonds_per_thread_target
                end_point -= 1
                n_bonds_per_tid -= n_family_members[end_point]
            end
        end
        balancedist[tid] = start_point:end_point
        start_point = end_point + 1
    end
    return balancedist
end

#
n_threads = 64
n_family_members = copy(sp.n_family_members)
point_ddist = BBVV.defaultdist(length(n_family_members), n_threads)
balanced_dist = balancedist(point_ddist, n_family_members)
workload_of_dist(balanced_dist, n_family_members)
