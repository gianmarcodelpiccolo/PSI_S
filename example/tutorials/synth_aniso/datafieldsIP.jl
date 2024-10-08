 
function build_obslist()
    IP_obs = Vector{obsinfo}()
    # Positional Inputs:
    # 1. Initialized observation vector
    # 2. Observation type (not used for teleseismic delay inversions)
    # 3. Relative path to data file
    # 4. Demean data?
    # 5. Solve for event statics?
    # 6. Solve for station statics?
    # 7. Initial noise guess for data
    # 8. Forward function for making predictions
    add_obs(IP_obs,"P_tts",string((@__DIR__),"/input/tts_P.dat"),false,false,false,0.125,"tt_P_Vp_fp_psi")

    return IP_obs
end

function build_fieldslist()
    IP_fields = Vector{fieldinfo}()
    # Positional Inputs:
    # 1. Initialized parameter fields vector
    # 2. Field name
    # 3. prior choice
    # 4. [min., max.] field values limits
    # 5. [min., max.] Latitude limits (deg.)
    # 6. [min., max.] Longitude limits (deg.)
    # 7. [min., max.] Elevation limits (km)
    # 8. [min., max.] Time limits for 4D
    # 9. reference value for non-random initialization, extrapolation value when squeezing is enabled
    add_field(IP_fields,"Vp","uniform",[5.8,7.2],[-0.3,0.3],[-0.3,0.3],[-0.01,0.01],[0.0,0.0],6.0)
    add_field(IP_fields,"fp","uniform",[0.0,0.12],[-0.3,0.3],[-0.3,0.3],[-0.01,0.01],[0.0,0.0],0.0)
    add_field(IP_fields,"psi","uniform",[-π,π],[-0.3,0.3],[-0.3,0.3],[-0.01,0.01],[0.0,0.0],0.0)

    return IP_fields 
end

