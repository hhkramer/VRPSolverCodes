function build_model(data::DataClsp)

   # T = tasks(data)
   # K = machines(data)
   # Q = capacity(data) # Q[k] is the capacity of the machine k

   I = items(data)
   T = periods(data)
   C = capacity(data)

   clsp = VrpModel()

   @variable(gap.formulation, x[t in T, k in K], Bin)
   @objective(gap.formulation, Min, sum(c(data,t,k) * x[t,k] for t in T, k in K))
   @constraint(gap.formulation, cov[t in T], sum(x[t,k] for k in K) == 1)

   # Build the directed graph Gᵏ=(Vᵏ,Aᵏ)
   function build_WagnerWhitin_graph(i::Int)
      max_inv = 0
      for t in 1:length(T)
         max_inv += demand(data, i, t)
      end

      V = Array{Int}(undef, length(T), length(T)+1)
      S = Array{Int}(undef, length(T), length(T)+1)

      fill!(S, -1)
      for period in 1:length(T)+1
         S[1, t] = 0
      end
      
      for period in 2:length(T)
         increment = 0
         max_inv_level = length(T) - period + 2
         for inv_level in 2:max_inv_level
            S[inv_level, period] = sum_demands(data, i, inv_level, inv_level + increment)
            increment += 1
         end
      end


      
   end

   function build_graph(i::Int)
      V = [v for v in 0:length(T)]

      v_source, v_sink = 0, length(T)
      L, U = 0, 1
      G = VrpGraph(gap, V, v_source, v_sink, (L,U))

      cap_res_id = add_resource(G, main = true) # R = Rᴹ = {cap_res_id}
      arc_plus_ids = [] # arc_plus_ids[t] is the id for the arc (t-1,t) with positive resource consumption

      # Build A
      for t in 1:length(T)
         arc_plus_id = add_arc(G, t-1, t) 
         add_arc_var_mapping(G, arc_plus_id, x[t,k])   
         set_arc_consumption(G, arc_plus_id, cap_res_id, w(data,t,k))
         push!(arc_plus_ids, [arc_plus_id]) # add arc id for the task t
 
         arc_minus_id = add_arc(G, t-1, t) 
         set_arc_consumption(G, arc_minus_id, cap_res_id, 0.0)         
      end

      # Accumulated resource consumption interval [l_v, u_v] for each vertex v \in V
      for v in V
         l_v, u_v = 0.0, Float64(Q[k])
         set_resource_bounds(G, v, cap_res_id, l_v, u_v)
      end
      return G, arc_plus_ids
   end

   graphs, packing_set_arcs = [], [[] for t in T]
   for k in K
      G, arc_plus_ids = build_graph(k) 
      add_graph(gap, G)
      push!(graphs, G)
      packing_set_arcs = [vcat(packing_set_arcs[t], arc_plus_ids[t]) for t in 1:length(T)]
   end

   set_arc_packing_sets(gap, [[(graphs[k],packing_set_arcs[t][k]) for k in K] for t in T]) 
   set_branching_priority(gap, "x", 1)
   
   return (gap, x)
end
