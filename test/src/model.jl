function build_WagnerWhitin_model(data::DataClsp, G::Array{InputGraph})
    I = items(data)
    T = periods(data)

    ww_model = Model(with_optimizer(Gurobi.Optimizer))

    @variable(ww_model, x[i in I, e in 1:length(G[i].A)], Bin)

    @objective(ww_model, Min, sum(G[i].A[e].cost * x[i, e] for i in I, e in 1:length(G[i].A)))

    @constraint(ww_model, source[i in I], sum(x[i, e] for e in leaving_arcs(G[i], G[i].V[1])) == 1)
    @constraint(ww_model, sink[i in I], sum(x[i, e] for e in entering_arcs(G[i], G[i].V[length(G[i].V)])) == 1)
    @constraint(ww_model, flow_cons[i in I, v in build_VV(data, G[i])], 
                sum(x[i, e] for e in entering_arcs(G[i], G[i].V[v])) - sum(x[i, e] for e in leaving_arcs(G[i], G[i].V[v])) == 0)


    print(ww_model)

    optimize!(ww_model)

    println(objective_value(ww_model))

    for i in I
        for e in 1:length(G[i].A)
            if value(x[i, e]) > 0.9
                println("x[$i, $e] = ", value(x[i, e]))
            end
        end
    end

end