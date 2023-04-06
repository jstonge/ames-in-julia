using Pkg; Pkg.activate("../../");
using ArgParse, Distributions, StatsBase, OrdinaryDiffEq, RecursiveArrayTools, DataFrames, SQLite

include("helpers.jl")

function parse_commandline()
  s = ArgParseSettings()

  @add_arg_table! s begin
      "--db"
      help = "Use Database to query parameters"
      "-L"
      arg_type = Int
      default = 5
      help = "LIMIT of rows"
      "-O"
      arg_type = Int
      default = 0
      help = "The OFFSET clause after LIMIT specifies how many rows to skip at the beginning of the result set."
      "--beta"
      arg_type = Float64
      default = 0.07
      help = "Spreading rate from non-adopter to adopter beta"
      "-a"
      arg_type = Float64
      default = 0.5
      help = "Negative benefits alpha"
      "-g"
      arg_type = Float64
      default = 1.
      help = "Recovery rate gamma, i.e. rate at which adopters loose behavioral trait"
      "-r"
      arg_type = Float64
      default = 0.1
      help = "Global behavioral diffusion ρ (allows the behaviour to spread between groups)"
      "-b"
      arg_type = Float64
      default = 0.18
      help = "Group benefits b"
      "-c"
      arg_type = Float64
      default = 1.05
      help = "Institutional cost c"
      "-m"
      arg_type = Float64
      default = 1e-4
      help = "Noise u"
      "-o"
      default = "."
      help = "Output file for results"
    end

  return parse_args(s)
end

function initialize_u0(;n::Int=20, L::Int=6, M::Int=20, p::Float64=0.01)
    G = zeros(L, n+1)
    for ℓ in 1:L
      for _ in 1:(floor(M/L))
        i = sum(rand(Binomial(1, p), n)) # how many total adopters?
        G[ℓ, i+1] += 1                   # everytime combination G[ℓ,i], count +1
      end
    end
  
    G = G ./ sum(G) # normalized by tot number of groups
    
    # ArrayPartition are nice because we can still access the level such as x[ℓ][i]
    return ArrayPartition(Tuple([G[ℓ,:] for ℓ=1:L]))
end

function source_sink!(du, u, p, t)
    G, L, n = u, length(u.x), size(u.x[1])[1]
    β, γ, ρ, b, c, μ = p
    Z, pop, R = zeros(L), zeros(L), 0.

    # Calculate mean-field coupling and observed fitness landscape
    for ℓ in 1:L
      n_adopt = collect(0:(n-1))
      Z[ℓ]    = sum(exp.(b*n_adopt .- c*(ℓ-1)) .* G.x[ℓ]) 
      pop[ℓ]  = sum(G.x[ℓ])
      R      += sum(ρ * n_adopt .* G.x[ℓ]) 
      pop[ℓ] > 0.0 && ( Z[ℓ] /= pop[ℓ] )
    end

    for ℓ = 1:L, i = 1:n
      n_adopt, gr_size = i-1, n-1

      # Diffusion events
      du.x[ℓ][i] = -γ*n_adopt*G.x[ℓ][i] - (ℓ-1)*β*(n_adopt+R)*(gr_size-n_adopt)*G.x[ℓ][i]
 
      n_adopt > 0 && ( du.x[ℓ][i] += β*(ℓ-1)*(n_adopt-1+R)*(gr_size-n_adopt+1)*G.x[ℓ][i-1])
      n_adopt < gr_size && ( du.x[ℓ][i] +=  γ*(n_adopt+1)*G.x[ℓ][i+1] )
 
      # Group selection process
      ℓ > 1 && ( du.x[ℓ][i] += ρ*G.x[ℓ-1][i]*(Z[ℓ] / Z[ℓ-1] + μ) - ρ*G.x[ℓ][i]*(Z[ℓ-1] / Z[ℓ]+μ) )
      ℓ < L && ( du.x[ℓ][i] += ρ*G.x[ℓ+1][i]*(Z[ℓ] / Z[ℓ+1] + μ) - ρ*G.x[ℓ][i]*(Z[ℓ+1] / Z[ℓ]+μ) )
    end
end

function run_source_sink(p)
  n, M = 20, 1000
  u₀ = initialize_u0(n=n, L=6, M=M, p=0.05)
  tspan = (1.0, 4000)
  
  # Solve problem
  prob = ODEProblem(source_sink!, u₀, tspan, p)
  return solve(prob, DP5(), saveat = 1., reltol=1e-8, abstol=1e-8)
end


function main()
    # β, γ, ρ, b, c, μ = 0.07, 1., 0.1, 0.18, 1.05, 0.0001
    args = parse_commandline()
    
    if isnothing(args["db"])
      β = args["beta"]
      γ = args["g"]
      ρ = args["r"]
      b = args["b"]
      c = args["c"]
      μ = args["m"]
      
      p = [β, γ, ρ, b, c, μ]
      sol = run_source_sink(p)
      write_sol2txt("$(args["o"])/sourcesink1_$(join(p, "_")).txt", sol) 

    else
      
      db = SQLite.DB(args["db"])
      con = DBInterface.execute(db, """SELECT * from sourcesink1 LIMIT $(args["O"]), $(args["L"])""") |> DataFrame
    
      for row in eachrow(con)
        β = row["beta"]
        γ = row["gamma"]
        ρ = row["rho"]
        b = row["b"]
        c = row["cost"]
        μ = row["mu"]
      
        p = [β,  γ, ρ, b, c, μ]
        sol = run_source_sink(p)    
        write_sol2txt("$(args["o"])/sourcesink1_$(join(p, "_")).txt", sol) 
      end
    end  
end

main()

# prototyping -------------------------------------------------------------------------------

using Plots

n, M = 20, 1000
u₀ = initialize_u0(n=n, L=6, M=M, p=0.001)
t_max = 1000
tspan = (0., t_max)

β, γ, ρ, b, c, μ = 0.22, 0.9, 0.1, 0.18, 1.85, 0.0001
p  = [β, γ, ρ, b, c, μ]
prob1 = ODEProblem(source_sink!, u₀, tspan, p)
sol1 = solve(prob1, DP5(), saveat=1, reltol=1e-8, abstol=1e-8)

inst_level, inst_level_prop = parse_sol(sol1)  # params: β, γ, ρ, b, c, μ, δ

plot_scatter(inst_level, inst_level_prop)
plot_scatter(inst_level, inst_level_prop, plot_prop=false)

