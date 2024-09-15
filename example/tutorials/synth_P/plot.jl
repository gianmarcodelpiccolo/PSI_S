using Plots
include(ENV["PSI_S"]*"/src/dependencies.jl")
mkpath("output/synth_P/figs/")

# -- loads the ensemble of models
models = load("output/synth_P/synth_P_inv.jld","MarkovChain")

# -- defines longitude and latitude ranges for plotting
φ, θ = collect(range(deg2rad.(-0.3),deg2rad.(0.3),length=100)), collect(range(deg2rad.(-0.3),deg2rad.(0.3),length=100))

# -- for every location, we populate the ensemble of values from the models and compute statistics (average, standard deviation)
Vp_mean, Vp_std = zeros(Float64,length(φ),length(θ)), zeros(Float64,length(φ),length(θ))
for i in eachindex(φ), j in eachindex(θ)
    lon, lat = φ[i], θ[j]
    x,y,z = @cartesian([lat],[lon],[6371.0])
    x,y,z = x[1],y[1],z[1]
    Vp_vals = Float64[]
    # -- loop over the sampled models
    for model in models
        voronoi_Vp = model.fields[1]        # -- first field of the model (Vp Voronoi diagram)
        ind = v_dist(x, y, z, voronoi_Vp)   # -- index of the closest Voronoi cell to location (x,y,z)
        push!(Vp_vals,voronoi_Vp.v[ind])    # -- pushes value of the closest Voronoi cell into the ensemble
    end
    Vp_mean[i,j] = mean(Vp_vals)            # -- compute average for the location
    Vp_std[i,j] = std(Vp_vals)              # -- compute standard deviation for the location
end

# -- plots
heatmap(rad2deg.(φ),rad2deg.(θ),Vp_mean',xlabel="longitude", ylabel="latitude", c=cgrad(:rainbow, rev=true),aspect_ratio=:equal, clim=(6,8))
savefig("output/synth_P/figs/Vp_mean.png")
heatmap(rad2deg.(φ),rad2deg.(θ),Vp_std',xlabel="longitude", ylabel="latitude", cmap=:magma,aspect_ratio=:equal)
savefig("output/synth_P/figs/Vp_std.png")
