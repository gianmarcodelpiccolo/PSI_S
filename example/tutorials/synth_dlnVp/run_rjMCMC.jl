
using Distributed
addprocs(parse(Int64, ARGS[2]))
@everywhere include(ENV["PSI_S"]*"/src/dependencies.jl")

const wrk_dir = @__DIR__
const name = string(ARGS[1])
const chain_id = parse(Int64,string(ARGS[3]))
run_RJMCMC(wrk_dir, name, chain_id)

