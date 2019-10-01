function build_model(data::DataClsp, G::Array{InputGraph})
    I = items(data)
    T = periods(data)

    model = Model(with_optimizer(CPLEX.Optimizer))

    @variable(model, x[i in I, e in 1:length(G[i].A)], Bin)

    @objective(model, Min, sum(G[i].A[e].cost * x[i, e] for i in I, e in 1:length(G[i].A)))

    @constraint(model, source[i in I], sum(x[i, e] for e in leaving_arcs(G[i], G[i].V[1])) == 1)
    @constraint(model, sink[i in I], sum(x[i, e] for e in entering_arcs(G[i], G[i].V[length(G[i].V)])) == 1)
    @constraint(model, flow_cons[i in I, v in build_VV(data, G[i])], 
                sum(x[i, e] for e in entering_arcs(G[i], G[i].V[v])) 
                - sum(x[i, e] for e in leaving_arcs(G[i], G[i].V[v])) == 0)
    @constraint(model, capacity[t in T], 
                sum(G[i].A[e].res_consumption * x[i, e] for i in I, e in arcs_leaving_period(G[i], t)) <= data.cap)


    #print(model)

    optimize!(model)

    println(objective_value(model))

    for i in I
        for e in 1:length(G[i].A)
            if value(x[i, e]) > 0.9
                println("x[$i, $e] = ", value(x[i, e]))
            end
        end
    end

end