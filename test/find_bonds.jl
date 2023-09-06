@testitem "find bonds" begin
    pc = PointCloud(1,1,1,0.5)
    δ = 0.51
    n_threads = 4
    owned_points = BBVV.defaultdist(pc.n_points, n_threads)
    bonds, n_family_members, init_dists = BBVV.find_bonds(pc, δ, owned_points, n_threads)
    @test bonds == [
        2
        3
        5
        1
        4
        6
        1
        4
        7
        2
        3
        8
        1
        6
        7
        2
        5
        8
        3
        5
        8
        4
        6
        7
    ]
    @test n_family_members == fill(3, 8)
    @test init_dists == fill(0.5, 24)
end
