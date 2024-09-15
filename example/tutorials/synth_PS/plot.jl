using Plots
include(ENV["PSI_S"]*"/src/dependencies.jl")
mkpath("output/synth_PS/figs/")

# -- loads the ensemble of models
models = load("output/synth_PS/synth_PS_inv.jld","MarkovChain")

# -- defines longitude and latitude ranges for plotting
φ, θ = collect(range(deg2rad.(-0.3),deg2rad.(0.3),length=100)), collect(range(deg2rad.(-0.3),deg2rad.(0.3),length=100))

# -- for every location, we populate the ensemble of values from the models and compute statistics (average, standard deviation)
Vp_mean, Vp_std = zeros(Float64,length(φ),length(θ)), zeros(Float64,length(φ),length(θ))
Vp2Vs_mean, Vp2Vs_std = zeros(Float64,length(φ),length(θ)), zeros(Float64,length(φ),length(θ))
Vs_mean, Vs_std = zeros(Float64,length(φ),length(θ)), zeros(Float64,length(φ),length(θ))
Vp_Vp2Vs_cor, Vp_Vs_cor = zeros(Float64,length(φ),length(θ)), zeros(Float64,length(φ),length(θ))
for i in eachindex(φ), j in eachindex(θ)
    lon, lat = φ[i], θ[j]
    x,y,z = @cartesian([lat],[lon],[6371.0])
    x,y,z = x[1],y[1],z[1]
    Vp_vals = Float64[]
    Vp2Vs_vals = Float64[]
    for model in models
        voronoi_Vp = model.fields[1]            # -- first field of the model (Vp Voronoi diagram)
        ind = v_dist(x, y, z, voronoi_Vp)       # -- index of the closest Voronoi cell to location (x,y,z)
        push!(Vp_vals,voronoi_Vp.v[ind])        # -- pushes value of the closest Voronoi cell into the ensemble
        voronoi_Vp2Vs = model.fields[2]         # -- second field of the model (Vp2Vs Voronoi diagram)
        ind = v_dist(x, y, z, voronoi_Vp2Vs)    # -- index of the closest Voronoi cell to location (x,y,z)
        push!(Vp2Vs_vals,voronoi_Vp2Vs.v[ind])  # -- pushes value of the closest Voronoi cell into the ensemble
    end
    Vs_vals = @. Vp_vals / Vp2Vs_vals           # -- given the sampled values of Vp and Vp2Vs, Vs is computed as Vp/Vp2Vs
    Vp_mean[i,j] = mean(Vp_vals)                # -- calculates Vp values average in the location
    Vp_std[i,j] = std(Vp_vals)                  # -- calculates Vp values standard deviation in the location
    Vp2Vs_mean[i,j] = mean(Vp2Vs_vals)          # -- same follows for Vp2Vs and Vs...
    Vp2Vs_std[i,j] = std(Vp2Vs_vals)
    Vs_mean[i,j] = mean(Vs_vals)
    Vs_std[i,j] = std(Vs_vals)
    Vp_Vp2Vs_cor[i,j] = cor(Vp_vals,Vp2Vs_vals) # -- calculates correlation between Vp and Vp2Vs values
    Vp_Vs_cor[i,j] = cor(Vp_vals,Vs_vals)       # -- calculates correlation between Vp and Vs values
end

# -- plots
heatmap(rad2deg.(φ),rad2deg.(θ),Vp_mean',xlabel="longitude", ylabel="latitude", c=cgrad(:rainbow, rev=true),aspect_ratio=:equal, clim=(6,8))
savefig("output/synth_PS/figs/Vp_mean.png")
heatmap(rad2deg.(φ),rad2deg.(θ),Vp_std',xlabel="longitude", ylabel="latitude", cmap=:magma,aspect_ratio=:equal)
savefig("output/synth_PS/figs/Vp_std.png")

heatmap(rad2deg.(φ),rad2deg.(θ),Vp2Vs_mean',xlabel="longitude", ylabel="latitude", c=cgrad(:roma, rev=true),aspect_ratio=:equal, clim=(1.6,2.2))
savefig("output/synth_PS/figs/Vp2Vs_mean.png")
heatmap(rad2deg.(φ),rad2deg.(θ),Vp2Vs_std',xlabel="longitude", ylabel="latitude", cmap=:magma,aspect_ratio=:equal)
savefig("output/synth_PS/figs/Vp2Vs_std.png")

heatmap(rad2deg.(φ),rad2deg.(θ),Vs_mean',xlabel="longitude", ylabel="latitude", c=cgrad(:rainbow, rev=true),aspect_ratio=:equal, clim=(2.5,5.0))
savefig("output/synth_PS/figs/Vs_mean.png")
heatmap(rad2deg.(φ),rad2deg.(θ),Vs_std',xlabel="longitude", ylabel="latitude", cmap=:magma,aspect_ratio=:equal)
savefig("output/synth_PS/figs/Vs_std.png")

heatmap(rad2deg.(φ),rad2deg.(θ),Vp_Vp2Vs_cor',xlabel="longitude", ylabel="latitude", c=cgrad(:RdBu, rev=false), aspect_ratio=:equal,clim=(-1,1))
savefig("output/synth_PS/figs/Vp_Vp2Vs_cor.png")
heatmap(rad2deg.(φ),rad2deg.(θ),Vp_Vs_cor',xlabel="longitude", ylabel="latitude", c=cgrad(:RdBu, rev=false), aspect_ratio=:equal,clim=(-1,1))
savefig("output/synth_PS/figs/Vp_Vs_cor.png")
