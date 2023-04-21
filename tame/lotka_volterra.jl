using Pkg; Pkg.activate("..");
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

function lotka_volterra!(du, u, p, t)
  N1, N2 = size(u)
  K = size(u, 1)
  α, γ, β  = p

  for n₁=1:N1, n₂=1:N2
      n_prey, n_pred = n₁-1, n₂-1 #distinguish indices from pop

      du[n₁,n₂] = 0 # we need to start from somewhere
      du[n₁,n₂] -= γ*n_pred*u[n₁,n₂] # predator death output
      n₁ < N1 && ( du[n₁,n₂] -= α*n_prey*(K-n_prey)*u[n₁,n₂]/K ) # prey birth ouput
      n₂ < N2 && ( du[n₁,n₂] -= β*n_prey*n_pred*u[n₁,n₂] ) # predation output
      n_prey > 0 && ( du[n₁,n₂] += α*(n_prey-1)*(K-n_prey+1)*u[n₁-1,n₂]/K ) # prey birth input
      n₂ < N2 && (du[n₁,n₂] += γ*( n_pred+1)*u[n₁,n₂+1] ) # pred death input 
      (n_pred > 0 && n₁ < N1) && ( du[n₁,n₂] += β*(n_prey+1)*(n_pred-1)*u[n₁+1,n₂-1] ) #predation input
      end
end


t_max = 500
tspan = (0., t_max)
nb_of_states = 25
α, γ, β = 0.9, 0.4, 0.04
u₀ = zeros(nb_of_states, nb_of_states)
u₀[11,7] = 1.

p  = [α, γ, β]
prob = ODEProblem(lotka_volterra!, u₀, tspan, p)
sol = solve(prob, Tsit5(), saveat=1)


# -------------------------------- plot simple ------------------------------- #


p1=plot_lv_solution_simple(sol, dims=1)
p2=plot_lv_solution_simple(sol, dims=2)
p3=heatmap(log10.(sol.u[2]))
plot(p1, p2, layout=(2,1))


# ---------------------------------- plot 2 ---------------------------------- #


# taking the ∑_{i}^{N} i/n * x_{i}  / ∑_{i}^{N} x_{i}
processing_sol1(x, n) = sum((collect(0:(n-1)) / n) .* x) / sum(x) 

function parse_sol(;dims=1)
  N1, N2 = size(first(sol))
  tmax = length(sol)-1
  inst_level = Dict()
  inst_level_prop = Dict()
  for i=1:(dims==1 ? N1 : N2)
    values = []
    values_prop = []
    for t=2:tmax
      # t,i,dims=10,1,2
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
  N = N2 = length(res)
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

plot_time(res_prey, res_prey_prop, "prey", plot_prop=true)
plot_time(res_prey, res_prey_prop, "prey", l=:bottomleft)

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
savefig("test_lv.png")