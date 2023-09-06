using BBVV
using Base.Threads

##--
RESPATH = joinpath(@__DIR__, "results")
simname = "run1"
l = 1.0
Δx = 1/20
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
set_top = findall(p -> p[2] ≥ l/2-Δx, eachcol(pc.position))
set_bottom = findall(p -> p[2] ≤ -l/2+Δx, eachcol(pc.position))
bc_bottom = VelocityBC(t -> -v0, set_bottom, 2)
bc_top = VelocityBC(t -> v0, set_top, 2)
bcs = [bc_bottom, bc_top]
path = joinpath(RESPATH, simname)
ispath(path) && rm(path; recursive=true, force=true)
!ispath(path) && mkpath(path) # create the path if it doe

##--
sp, gs, vtls = BBVV.init_simulation(pc, mat, bcs, n_timesteps, export_freq, export_path)


##--
Δt = BBVV.calc_stable_timestep(sp, vtls)
Δt½ = 0.5Δt

##--
t = 1
time = t * Δt

##-- update velocity_half
@threads :static for i in 1:sp.pc.n_points
    for d in 1:3
        gs.velocity_half[d, i] = gs.velocity[d, i] + gs.acceleration[d, i] * Δt½
    end
end

##-- apply_bcs
@sync for bc in sp.bcs
    Threads.@spawn begin
        value = bc.fun(t)
        dim = bc.dim
        @inbounds @simd for i in bc.point_id_set
            gs.velocity_half[dim, i] = value
        end
    end
end

##-- update_disp_and_position
@threads :static for i in 1:sp.pc.n_points
    for d in 1:3
        gs.displacement[d, i] += gs.velocity_half[d, i] * Δt
        gs.position[d, i] += gs.velocity_half[d, i] * Δt
    end
end

##-- compute_forcedensity
@threads :static for tls in vtls

    for i in tls.tl_points
        for cbond in sp.hood_range[i]
            j = sp.bonds[cbond]
            L = sp.init_dists[cbond]
            ϑx = gs.position[1, j] - gs.position[1, i]
            ϑy = gs.position[2, j] - gs.position[2, i]
            ϑz = gs.position[3, j] - gs.position[3, i]
            l = sqrt(ϑx * ϑx + ϑy * ϑy + ϑz * ϑz)
            ε = (l - L) / L
            temp = gs.bond_failure[cbond] * ε / l
            tempi = temp * sp.mat.bc * sp.pc.volume[j]
            tempj = temp * sp.mat.bc * sp.pc.volume[i]
            tls.b_int[1, i] += tempi * ϑx
            tls.b_int[2, i] += tempi * ϑy
            tls.b_int[3, i] += tempi * ϑz

            # failure mechanism
            if ε > sp.mat.εc
                gs.bond_failure[cbond] = 0
            end
            gs.n_active_family_members[i] += gs.bond_failure[cbond]
        end
    end
end

##-- reduction
gs.b_int .= 0
gs.n_active_family_members .= 0
for tls in vtls
    for (gi, li) in tls.pointmap
        for d in 1:3
            gs.b_int[d, gi] += tls.b_int[d, li]
        end
        gs.n_active_family_members[gi] += tls.n_active_family_members[li]
    end
end

##-- calc_damage
@threads :static for i in 1:sp.pc.n_points
    gs.damage[i] = 1 - gs.n_active_family_members[i] / sp.n_family_members[i]
end

##-- compute equation of motion
@threads :static for i in 1:sp.pc.n_points
    for d in 1:3
        gs.acceleration[d, i] = (gs.b_int[d, i]) / sp.mat.rho
        gs.velocity[d, i] = gs.velocity_half[d, i] + gs.acceleration[d, i] * Δt½
    end
end

##-- export results
if mod(t, sp.export_freq) == 0
    export_vtk(sp, gs, t, time)
end
