using Pkg; Pkg.activate("../..");
using OrdinaryDiffEq, Plots, RecursiveArrayTools
using OrdinaryDiffEq:ODESolution
using Plots.PlotMeasures

function initialize_u0(; N::Int=25)
  G = zeros(N, N)
  G[12,8] = 1
  return G
end

function plot_lv_solution_simple(s; dims=2)
  t0=2
  ys = sum(sol.u[t0], dims=dims) # sum over rows (n2 axis) -> 25x1 matrix
  ys = dims == 1 ? transpose(ys) : ys
  xs = 1:length(ys)
  p = plot(xs, ys, marker = 2, label="t=2", legend = :outertopright)
  for t=5:3:34
    ys = sum(s.u[t], dims=dims) # sum over rows (n2 axis) -> 25x1 matrix
    ys = dims == 1 ? transpose(ys) : ys
    plot!(xs, ys, marker = 2, label="t=$(t)")
  end
  # type_ent = dims == 1 ? "non-programmers" : "programmers"
  type_ent = dims == 1 ? "prey" : "predator"
  xlabel!("Number of $(type_ent)")
  ylabel!("Occupation number")
  p    
end


function life_cycle_research_groups!(du, u, p, t)
    """
    Guillaume's model
    =================

     - s_m(t): prob for a node of membership m to be suceptible at time t
     - c_{n,i}(t) prob to observe i inf nodes within a group of size n at time t
    
     """
    C, S, N, I, M = u[1] u[2], size(u[1]), size(u[2], 1)
    μ, β = p
    r, ρ = zeros(I), zeros(I)
   
    # mean-field r(t)
    # p_n is the group size disribution??
    # n*p_n is the size n of a group drawn from p_n
    # for n=1:N, i=1:I
    # n_susceptible, n_inf = n-1, i-1
    # num_r += sum(β*i*(n_susceptible-n_inf)*C[n,i]*p[n])
    # denum_r += sum((n_susceptible-n_inf)*C[n,i]*p[n])
    # end
    # r = num_r / denum_r
  
    # mean-field ρ(t)
    # !todo What are we doing here? 
    # !todo Is it the equivalent of `Z[ℓ] = sum(exp.(b*[0,1,...n_adopt]) .- c*(ℓ-1)) .* G.x[ℓ])` in sourcesink_m2
    # for m=1:M
    # num_rho + =sum((n_membership-1)*n_membership*S[m]*g[m])
    # denum_rho += sum(n_membership*S[m]*g[m])
    # end
    # ρ = r * (num_rho / denum_rho)
    
    # for m=1:M
    #   n_membership = m - 1
    #   dS[m] = μ*(1 - S[m]) - n_membership*r*S[m]
    # end

    for n=1:N, i=1:I
      n_susceptible, n_inf = n-1, i-1   # we distinguish indices from actual values.
      
      dC[n,i] = 0
      
      dC[n,i] -= μ*n_inf*C[n,i]
      dC[n,i] -= (n_susceptible-n_inf)*(β*i + ρ)*C[n,i]

      n_inf < I && ( dC[n,i] += μ*(n_inf+1)*C[n,i+1] )
      n_inf > 0 && ( dC[n,i] += μ*(n_susceptible-n_inf+1)*(β*(n_inf-1))*C[n,i-1] )
    end
end


# μ  = 0.001   # inflow new students-non coders
# νₙ = 0.01   # death rate non-coders
# νₚ = 0.01    # death rate coders
# β  = 0.01     # conversion rates non-prog to prog
# p  = [μ, νₙ, νₚ, β]

# t_max = 100
# tspan = (0., t_max)
# nb_of_states = 25
# u₀ = initialize_u0(N=nb_of_states)

# prob = ODEProblem(life_cycle_research_groups!, u₀, tspan, p)
# sol = solve(prob, Tsit5(), saveat=1)

t_max = 500
tspan = (0., t_max)
nb_of_states = 25
α, γ, β = 0.9, 0.4, 0.04
u₀ = zeros(nb_of_states, nb_of_states)
u₀[11,7] = 1.

p  = [α, γ, β]
prob = ODEProblem(lotka_volterra!, u₀, tspan, p)
sol = solve(prob, Tsit5(), saveat=1)

# param_str = join(["$(x)$(y)" for (x,y) in zip(p, ["μ", "νₙ", "νₚ", "β"])], "_")


p1=plot_lv_solution_simple(sol, dims=1)
p2=plot_lv_solution_simple(sol, dims=2)
p3=heatmap(log10.(sol.u[2]))
plot(p1, p2, layout=(2,1))

processing_sol1(x, n) = sum((collect(0:(n-1)) / n) .* x) / sum(x) 

function parse_sol(;dims=1)
  N1, N2 = size(sol[1])
  tmax = length(sol)-1
  inst_level = Dict()
  inst_level_prop = Dict()
  for i=1:(dims==1 ? N1 : N2)
    values = []
    values_prop = []
    for t=2:tmax
      n = dims==1 ? length(sol[t][:,i]) : length(sol[t][i,:])
      x = dims==1 ? sol[t][:,i] : sol[t][i,:]
      out = processing_sol1(x,n)
      push!(values, out)
      out = sum(x)
      push!(values_prop, out)
    end
    inst_level[i] = values
    inst_level_prop[i] = values_prop
  end
  return inst_level, inst_level_prop
end

res_pred, res_pred_prop = parse_sol(dims=1)
res_prey, res_prey_prop = parse_sol(dims=2)


function plot_time(res, res_prop, type_ent; plot_prop=false, l=:topleft)
  N = length(res)
  pal = palette(:viridis, N)
  tmax = length(res[1])
  xs = 1:tmax
  
  ys = plot_prop ? [res_prop[i] for i in 1:N2] : [res[i] for i in 1:N2]
  p4 = plot(xs, ys[1], xaxis=:log10, c = pal[1], label="1", legendtitle= "# $(type_ent)", left_margin = 10mm, bottom_margin = 10mm)
  for n2 in 2:(N2-1)
    if Int(n2) ∈ (7, 15)
      plot!(xs, ys[Int(n2)], xaxis=:log10, c = pal[Int(n2)], label="$(n2)", legend=false)
    end 
    plot!(xs, ys[Int(n2)], xaxis=:log10, c = pal[Int(n2)], label="", legend=false)
  end
  p4 = plot!(xs, ys[N2], xaxis=:log10, c = pal[N2], label=N2, legendtitle= "# $(type_ent)", legend = l)
  
  if plot_prop == false
    global_freq = [sum([res[i][t]*res_prop[i][t] for i in 1:N]) for t in 1:tmax]
    plot!(1:tmax, global_freq, linestyle=:dash, color=:black, width = 1.5, label = "global") 
  end

  ylabel!(plot_prop ? "proportion $(type_ent)" : "frequency $(type_ent)")
  xlabel!("time")
  return p4  
end

p_tot = plot(
  plot_time(res_prey, res_prey_prop, "prey", l=:bottomleft),
  plot_time(res_prey, res_prey_prop, "prey", plot_prop=true),
  p1,
  plot_time(res_pred, res_pred_prop, "pred"),
  plot_time(res_pred, res_pred_prop, "pred", plot_prop=true),
  p2,
  layout= (2,3)
)
plot!(size=(1600,850))
savefig("prey_pred_quadrant.pdf")

plot(p1,p2,p3,p4)
plot!(size=(1200,700))
p1
savefig("test_lv.png")


ys = [res_prop[i] for i in 1:N2]
p4 = plot(xs, ys[1], xaxis=:log10, c = pal[1], label="1", legendtitle= "# prey")
for n2 in 2:(N2-1)
  plot!(xs, ys[Int(n2)], xaxis=:log10, c = pal[Int(n2)], label="", legend=false)
end
p4 = plot!(xs, ys[N2], xaxis=:log10, c = pal[N2], label=N2, legendtitle= "# prey", legend=true)
ylabel!("prop prey")
xlabel!("time")



# -------------------- Parse sol global average group size ------------------- #


N1, N2 = size(sol[1])
maxsize = (N1-1)+(N2-1)
I = zeros(maxsize, length(sol.t))
weighted_avg = zeros(maxsize, length(sol.t))
size_dis = zeros(maxsize, length(sol.t))
dist_gsize = zeros(maxsize, length(sol.t))
sommeGNI = zeros(maxsize, length(sol.t))

for t=1:length(sol.t)
  # t=20
  for n1=1:N1
    for n2=1:N2
      # t,n1,n2 = 2,1,3
      coder, noncoder = n1-1, n2-1 
      gsize = coder+noncoder
      
      G_nil = sol[t][n1,:] # sol of group with n non-coders at time t
      gsize > 0 && ( dist_gsize[gsize, t] += G_nil[n2] )  
      gsize > 0 && ( weighted_avg[gsize, t] += (coder/gsize)*G_nil[n2] )
      gsize > 0 && ( size_dis[gsize, t] += G_nil[n2] )
    end
  end
  
  for gsize=1:maxsize
      I[gsize,t] = weighted_avg[gsize, t]/size_dis[gsize, t]
  end
  
  for i=1:t_max
    dist_gsize[:,i] = dist_gsize[:,i] / sum(dist_gsize[:,i])
  end

end

dist_gsize

p = scatter(dist_gsize[1,1:10], legendtitle="grsize", legend=:outertopright, label="1", scale=:log10)
for i=2:nb_of_states
  scatter!(dist_gsize[i,10:t_max], label="$(i)", scale=:log10)
end
p

