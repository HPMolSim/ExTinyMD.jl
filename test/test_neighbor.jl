using CellListMap
using StaticArrays


T = Float32;
x = [SVector(10 * rand(T), 10 * rand(T), 10 * rand(T)) for _=1:100];
@time celllist = InPlaceNeighborList(x = x, cutoff = T(1.0), unitcell = [T(10.0), T(10.0), T(10.0)], parallel = false);
@time nlist = neighborlist!(celllist);

x_new = [SVector(10 * rand(Float64), 10 * rand(Float64), 10 * rand(Float64)) for _=1:100];
@time update!(celllist, x_new);
@time neighborlist!(celllist);