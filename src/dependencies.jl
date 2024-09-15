# Package Dependencies
using JLD
using Interpolations
using DelimitedFiles
using IntervalSets
using Statistics
using LazySets
using Random
using DataStructures
using NearestNeighbors
using TauP
using Dierckx
using Distributed
using DistributedArrays
using SharedArrays
using MiniQhull
using LinearAlgebra
using SparseArrays

# Global Constants
const R = 6371.0

# Package Files
include("geodesics.jl")
include("structures.jl")
include("initializeIP.jl")
include("voronoi_utilities.jl")
include("RJMCMC.jl")
include("statics_solver.jl")
include("parallel_tempering.jl")
include("delayed_rejection.jl")
include("forward_functions/utilities.jl")
include("forward_functions/P_phases/forward.jl")
include("forward_functions/S_phases/forward.jl")
# This ray tracer should really be a stand-alone package
include("raytracer/grid.jl")
include("raytracer/dijkstra.jl")
include("raytracer/bfm.jl")
include("raytracer/shortest_path.jl")
include("raytracer/chull_carving.jl")
include("raytracer/relocations.jl")
