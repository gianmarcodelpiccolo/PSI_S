using Plots
using DirectionalStatistics
using IntervalSets

include(ENV["PSI_S"]*"/src/dependencies.jl")
mkpath("output/synth_aniso/figs/")

function arrow0!(x, y, u, v; s=1, lc=:black, la=1)
    u,v = u/s, v/s
    plot!([x-u/2,x+u/2], [y-v/2,y+v/2], lc=lc, la=la, legend=false)
end

function updateT(T,f,𝛙,γ,N)
    i = pi/2 - γ
    ϕ = 3*pi/2 - 𝛙

    T[1,1] += f*sin(i)^2*sin(ϕ)^2/N
    T[1,2] += f*sin(i)^2*sin(ϕ)*cos(ϕ)/N
    T[1,3] += f*sin(i)*cos(i)*sin(ϕ)/N

    T[2,1] += f*sin(i)^2*sin(ϕ)*cos(ϕ)/N
    T[2,2] += f*sin(i)^2*cos(ϕ)^2/N
    T[2,3] += f*sin(i)*cos(i)*cos(ϕ)/N

    T[3,1] += f*sin(i)*cos(i)*sin(ϕ)/N
    T[3,2] += f*sin(i)*cos(i)*cos(ϕ)/N
    T[3,3] += f*cos(i)^2/N
end

# -- loads the ensemble of models
models = load("output/synth_aniso/synth_aniso_inv.jld","MarkovChain")

# -- defines longitude and latitude ranges for plotting
φ, θ = collect(range(deg2rad.(-0.3),deg2rad.(0.3),length=100)), collect(range(deg2rad.(-0.3),deg2rad.(0.3),length=100))

# -- for every location, we populate the ensemble of values from the models and compute statistics (average, standard deviation)
Vp_mean, Vp_std = zeros(Float64,length(φ),length(θ)), zeros(Float64,length(φ),length(θ))
fp_mean, fp_std = zeros(Float64,length(φ),length(θ)), zeros(Float64,length(φ),length(θ))
psi_mean, psi_std = zeros(Float64,length(φ),length(θ)), zeros(Float64,length(φ),length(θ))
v1,v2 = zeros(Float64,length(φ),length(θ)), zeros(Float64,length(φ),length(θ)) 
Vp_fp_cor = zeros(Float64,length(φ),length(θ))
v1_T, v2_T = zeros(Float64,length(φ),length(θ)), zeros(Float64,length(φ),length(θ)) 
fp_T = zeros(Float64,length(φ),length(θ))
DB = zeros(Float64,length(φ),length(θ))
λ1, λ2 = zeros(Float64,length(φ),length(θ)), zeros(Float64,length(φ),length(θ))
T = zeros(Float64,3,3)
for i in eachindex(φ), j in eachindex(θ)
    lon, lat = φ[i], θ[j]
    x,y,z = @cartesian([lat],[lon],[6371.0])
    x,y,z = x[1],y[1],z[1]
    Vp_vals = Float64[]
    fp_vals = Float64[]
    psi_vals = Float64[]
    T .= 0.0
    for model in models
        voronoi_Vp = model.fields[1]       # -- first field of the model (Voronoi diagram)
        ind = v_dist(x, y, z, voronoi_Vp)  # -- index of the closest Voronoi cell
        push!(Vp_vals,voronoi_Vp.v[ind])
        voronoi_fp = model.fields[2]      
        ind = v_dist(x, y, z, voronoi_fp)  
        push!(fp_vals,voronoi_fp.v[ind])
        voronoi_psi = model.fields[3]      
        ind = v_dist(x, y, z, voronoi_psi) 
        push!(psi_vals,voronoi_psi.v[ind])
        updateT(T,fp_vals[end],psi_vals[end],0.0,length(models))
    end
    Vp_mean[i,j], Vp_std[i,j] = mean(Vp_vals), std(Vp_vals)
    fp_mean[i,j], fp_std[i,j] = mean(fp_vals), std(fp_vals)
    psi_mean[i,j], psi_std[i,j] = Circular.mean(psi_vals, 0..π), Circular.std(psi_vals, 0..π) # -- !!!! circular statistics !!!!
    v1[i,j] = fp_mean[i,j] * cos(psi_mean[i,j]) * 1
    v2[i,j] = fp_mean[i,j] * sin(psi_mean[i,j]) * 1
    Vp_fp_cor[i,j] = cor(Vp_vals,fp_vals)
    vecs = eigvecs(T)
    vals = abs.(eigvals(T))
    inds = sortperm(vals)
    ind_a, ind_b, ind_c = reverse(inds)
    HexAxis = [-(vals[ind_a])*vecs[3,ind_a],(vals[ind_a])*vecs[1,ind_a],(vals[ind_a])*vecs[2,ind_a]]
    v1_T[i,j], v2_T[i,j] = HexAxis[2], HexAxis[3]
    fp_T[i,j] = sqrt(v1_T[i,j]^2+v2_T[i,j]^2)
    λ1[i,j], λ2[i,j] = abs(vals[ind_a]), abs(vals[ind_b])
    DB[i,j] = sqrt(0.5)*sqrt(((λ1[i,j]^2-λ2[i,j]^2))/(λ1[i,j]^2+λ2[i,j]^2))
end

# -- plots
Plots.heatmap(rad2deg.(φ),rad2deg.(θ),Vp_mean',xlabel="longitude", ylabel="latitude", c=cgrad(:rainbow, rev=true), aspect_ratio=:equal, clim=(5.5,7.5)) 
Plots.savefig("output/synth_aniso/figs/Vp_mean.png")
Plots.heatmap(rad2deg.(φ),rad2deg.(θ),Vp_std',xlabel="longitude", ylabel="latitude", cmap=:magma, aspect_ratio=:equal)
Plots.savefig("output/synth_aniso/figs/Vp_std.png")

φp, θp = Float64[], Float64[]
v1p, v2p = Float64[], Float64[]
for i in eachindex(φ)[1:5:end], j in eachindex(θ)[1:5:end]
    push!(φp,φ[i])
    push!(θp,θ[j])
    push!(v1p,v1[i,j])
    push!(v2p,v2[i,j])
end
scale = 8
Plots.heatmap(rad2deg.(φ),rad2deg.(θ),fp_mean',xlabel="longitude", ylabel="latitude", c=cgrad(:davos), aspect_ratio=:equal,  clim=(0,0.1))
for i in eachindex(φp)
    arrow0!(rad2deg.(φp[i]), rad2deg.(θp[i]), v1p[i], v2p[i]; s=7, lc=:black, la=1)
end
Plots.savefig("output/synth_aniso/figs/fp_anisotropy.png")

Plots.heatmap(rad2deg.(φ),rad2deg.(θ),fp_std',xlabel="longitude", ylabel="latitude", cmap=:magma, aspect_ratio=:equal)
Plots.savefig("output/synth_aniso/figs/fp_std.png")

Plots.heatmap(rad2deg.(φ),rad2deg.(θ),rad2deg.(psi_std)',xlabel="longitude", ylabel="latitude", cmap=:magma, aspect_ratio=:equal)
Plots.savefig("output/synth_aniso/figs/psi_std.png")

Plots.heatmap(rad2deg.(φ),rad2deg.(θ),Vp_fp_cor',xlabel="longitude", ylabel="latitude", c=cgrad(:RdBu, rev=false), aspect_ratio=:equal,clim=(-1,1))
Plots.savefig("output/synth_aniso/figs/Vp_fp_cor.png")

φp, θp = Float64[], Float64[]
v1p, v2p = Float64[], Float64[]
λ1p, λ2p = Float64[], Float64[]
for i in eachindex(φ)[1:6:end], j in eachindex(θ)[1:6:end]
    push!(φp,φ[i])
    push!(θp,θ[j])
    push!(v1p,v1_T[i,j])
    push!(v2p,v2_T[i,j])
    push!(λ1p,λ1[i,j])
    push!(λ2p,λ2[i,j])
end
t = range(0,2π,length=100)
Plots.heatmap(rad2deg.(φ),rad2deg.(θ),fp_T',xlabel="longitude", ylabel="latitude", c=cgrad(:davos), aspect_ratio=:equal,  clim=(0,0.1))
for i in eachindex(φp)
    # arrow0!(rad2deg.(φp[i]), rad2deg.(θp[i]), v1p[i], v2p[i]; s=scale, lc=:black, la=1)
    ψ = atan(v2p[i],v1p[i])
    x = @. (sqrt(λ1p[i])*cos(ψ)*cos(t) - sqrt(λ2p[i])*sin(ψ)*sin(t))/scale + rad2deg.(φp[i])
    y = @. (sqrt(λ1p[i])*sin(ψ)*cos(t) + sqrt(λ2p[i])*cos(ψ)*sin(t))/scale + rad2deg.(θp[i])
    Plots.plot!(x,y,c=:black,legend=false)
end
Plots.savefig("output/synth_aniso/figs/fp_anisotropy_tensor.png")

Plots.heatmap(rad2deg.(φ),rad2deg.(θ),DB',xlabel="longitude", ylabel="latitude", c=cgrad(:magma, rev=true), aspect_ratio=:equal,  clim=(0,1))
Plots.savefig("output/synth_aniso/figs/DB.png")







