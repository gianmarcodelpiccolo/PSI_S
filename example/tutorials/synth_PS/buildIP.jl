# Required Environment Variables
ENV["TAUP_JAR"] = ENV["PSI_S"]*"/src/TauP_Toolkit/TauP-2.7.0-SNAPSHOT5/lib/TauP-2.7.0-SNAPSHOT5.jar"
# Local Inversion Directory
wrk_dir = @__DIR__
include(ENV["PSI_S"]*"/src/dependencies.jl")
include("datafieldsIP.jl")

################################
## -- parameters selection -- ##
################################

# -- DEFINE GENERAL PARAMETERS FOR SEISMIC IMAGING
IP = IPConst(
    "synth_PS",                         # -- name of the run
	"input/evt.dat",           # -- events data file 
	"input/sta.dat",           # -- stations data file 
	"input/homogeneous.dat",   # -- 1D velocity model                            
	DomainGeoBoundaries(        # -- GEOGRAPHIC BOUNDARIES FOR THE INVERSION DOMAIN
        [-0.3,0.3],             # -- latitude limits [°]
        [-0.3,0.3],             # -- longitude limits [°]
        [0.0,0.0]               # -- depth limits [km]
    ),      
    MonteCarloSolver(           # -- PARAMETERS FOR RJMCMC SOLVER
        2,                      # -- misfit function L-n norm (L1, L2...)
        4e4,                    # -- total number of iterations
        12,                     # -- number of gianMarkov chains
        "linear",               # -- Nearest-Neighbour(NN) interpolation algorithm ("linear" -> linear search,"kdtree" -> NearestNeighbors.jl [deprecated])
        0.1,                    # -- perturb.size (fraction of lims provided for each voronoi diagram in datafieldsIP.jl)
        0,                      # -- initial step-1 range
        "uniform",              # -- prior choice for the number of nuclei ("uniform","log-uniform")   
        10000,                  # -- maximum number of nuclei [memory usage depends on this limit]
        100,                    # -- initial number of nuclei
        true,                   # -- hierarchical Bayesian inversion (include noise random variables in posterior definition for each obs.)
        true,                   # -- initializes randomly models in chain
        false,                  # -- squeezing (outside Voronoi diagrams the reference values provided in datafieldsIP.jl are used)
        1e2,                    # -- interval to print out on progress' state 
        1e3                     # -- saving interval [chain thinning] (only one model every n iterations is stored)
    ),   
    DelayedRejection(           # -- PARAMETERS FOR DELAYED REJECTION ALGORITHM  (only for rjMcMC steps 1 and 2)
        false,                  # -- status
        4                       # -- proposal st.dev rescaling
    ),
    ParallelTempering(          # -- PARAMETERS FOR PARALLEL TEMPERING
        false,                  # -- status
        1e3,                    # -- first P.T. swaps after n iterations
        1.0,                    # -- cold chains temperature
        1,                      # -- number of cold chains (T=1)
        10,                     # -- hottest chain temperature
        1e3,                    # -- swaps every n iterations
        true                    # -- increase prosal perturbation size as sqrt(T) 
    ),  
    StaticCorrections(          # -- PARAMETERS FOR EVENT/STATION STATIC CORRECTIONS INVERSION (LSQR algorithm)
        1e4,                    # -- first static evaluation after n iteration (set 1 to run at the beginning)
        1e4,                    # -- static evaluatins runs every n iteration
        1e6                     # -- expected ratio station-statics/residual uncertainties (damping parameter for station statics)
    ),
    RayTracing(                 # -- PARAMETERS FOR RAY-TRACING (TELESEISMIC AND LOCAL EARTHQUAKES DATA)
        0.4,                    # -- ray-paths discretization length [km]
        "linear",               # -- interpolation along ray-paths ["linear","spline"]
        "iasp91",               # -- TauP velocity model used for teleseismic ray-tracing (if used)
        false,                  # -- forces ray-tracing using TauP for all the rays
        6380.0,                 # -- maximum Earth's radius according to TauP extended velocity model
        [300,300,1],            # -- shortest-path grid nodes along each coordinate (latitude, longitude, radius)
        5,                      # -- shortest-path grid forward star level
        false,                  # -- allows noise perturbation to grid nodes 
        false,                  # -- grid-carving, might reduce significantly shortest-path execution time
        1e4                     # -- ray-tracing running every n iterations
    ),                                        
    EarthquakesRelocation(      # -- PARAMETERS FOR EARTHQUAKES RELOCATION (ONLY FOR 3D LET) [EXPERIMENTAL]
        false,                  # -- status
        3                       # -- relocation runs every n ray-tracing
    ),
    Bayesian4Dimaging(          # -- PARAMETERS FOR 4D IMAGING [EXPERIMENTAL]
        false                   # -- status
    )
)     


##########################
## -- initialization -- ##
##########################

print("\ninitializing the inverse problem...\n")
IP_fields = build_fieldslist()
IP_obs = build_obslist()
rnodes, rays, evtsta, observables, LocalRaysManager = initialize_IP(IP, IP_obs, IP_fields)

mkpath("output/"*IP.name)
save(
    string("output/", IP.name, "/", IP.name, ".jld"), "IP", IP,
    "IP_fields", IP_fields, 
    "rnodes", rnodes, 
    "rays", rays, 
    "evtsta", evtsta, 
    "observables", observables,
    "LocalRaysManager", LocalRaysManager
)

# Save copy of input parameters
# cp("buildIP.jl", "output/"*IP.name*"/buildIP.jl")
# cp("datafieldsIP.jl", "output/"*IP.name*"/datafieldsIP.jl")


