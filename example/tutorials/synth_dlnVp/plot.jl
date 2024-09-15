using Plots
include(ENV["PSI_S"]*"/src/dependencies.jl")
mkpath("output/synth_dlnVp/figs/")

# -- loads the ensemble of models
models = load("output/synth_dlnVp/synth_dlnVp_inv.jld","MarkovChain")

# -- defines longitude and latitude ranges for plotting
φ, θ = collect(range(deg2rad.(-0.3),deg2rad.(0.3),length=100)), collect(range(deg2rad.(-0.3),deg2rad.(0.3),length=100))

# -- for every location, we populate the ensemble of values from the models and compute statistics (average, standard deviation)
dlnVp_mean, dlnVp_std = zeros(Float64,length(φ),length(θ)), zeros(Float64,length(φ),length(θ))
for i in eachindex(φ), j in eachindex(θ)
    lon, lat = φ[i], θ[j]
    x,y,z = @cartesian([lat],[lon],[6371.0])
    x,y,z = x[1],y[1],z[1]
    dlnVp_vals = Float64[]
    for model in models
        voronoi_dlnVp = model.fields[1]           # -- first field of the model (Vp Voronoi diagram)
        ind = v_dist(x, y, z, voronoi_dlnVp)      # -- index of the closest Voronoi cell to location (x,y,z)
        push!(dlnVp_vals,voronoi_dlnVp.v[ind])    # -- pushes value of the closest Voronoi cell into the ensemble
    end
    dlnVp_mean[i,j] = mean(dlnVp_vals)      # -- compute average for the location
    dlnVp_std[i,j] = std(dlnVp_vals)        # -- compute standard deviation for the location
end

# -- plots
heatmap(rad2deg.(φ),rad2deg.(θ),dlnVp_mean',xlabel="longitude", ylabel="latitude", c=cgrad(:seismic, rev=true),aspect_ratio=:equal,clims=(-0.3,0.3))
savefig("output/synth_dlnVp/figs/dlnVp_mean.png")
heatmap(rad2deg.(φ),rad2deg.(θ),dlnVp_std',xlabel="longitude", ylabel="latitude", cmap=:magma,aspect_ratio=:equal)
savefig("output/synth_dlnVp/figs/dlnVp_std.png")
