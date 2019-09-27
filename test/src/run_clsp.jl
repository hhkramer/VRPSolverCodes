# using VrpSolver, JuMP
using JuMP, Gurobi

include("data_clsp.jl")
# include("model_clsp.jl")
# include("param_clsp.jl")
include("model.jl")

appfolder = dirname(@__FILE__)

# if length(ARGS) == 0
   #  cutoff = 12o44
    inputfile = "$appfolder/../data/X11119A"#testSmall"
# else
#     inputfile = ARGS[1] #reading the first parameter passed through command line
#     if (length(ARGS) == 1)
#       cutoff = 100000.0
#    else
#       cutoff = parse(Float64,ARGS[2]) #reading the second parameter passed through command line
#    end
# end

data = read_clsp_data(inputfile)

@show data

G = InputGraph[]
for i in 1:data.numItems
   # V = build_V(data, i)

   # for v in V
   #    println(v.vertex_id, " ", v.i, " ", v.s, " ", v.t, " ", v.inv_qtd)
   # end

   # A = build_A(data, i, V)

   # for a in A
   #    println(a.arc_id, " ", a.i, " ", a.origin, " ", a.destination, " ", a.prod_qtd, " ", a.res_consumption, " ", a.cost)
   # end
   g = build_G(data, i)
   push!(G, g)
   println(sizeof(G))
end

model = build_WagnerWhitin_model(data, G)


# (model, x) = build_model(data)
# optimizer = VrpOptimizer(model, "$appfolder/../config/GAP.cfg")
# set_cutoff(optimizer, cutoff)

# status = optimize(optimizer)

# Output
# println("Status is $status")
# if status == :Optimal
#    println("Objective value : $(get_objective_value(optimizer))")
#    for t in data.tasks
#       for k in data.machines
#          println("x[$t,$k] : $(get_value(optimizer, x[t,k]))")
#       end
#    end
# end
