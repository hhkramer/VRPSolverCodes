function build_WagnerWhitin_model(data::DataClsp, G::Array{InputGraph})
    I = items(data)
    T = periods(data)

    ww_model = Model(with_optimizer(Gurobi.Optimizer))

    @variable(ww_model, x[i in I, e in 1:length(G[i].A)], Bin)

    @objective(ww_model, Min, sum(G[i].A[e].cost * x[i, e] for i in I, e in 1:length(G[i].A)))

    @constraint(ww_model, source[i in I], sum(x[i, e] for e in leaving_arcs(G[i], G[i].V[1])) == 1)


    print(ww_model)

end