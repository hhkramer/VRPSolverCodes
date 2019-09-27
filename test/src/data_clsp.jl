import Base.show, Base.print

using DelimitedFiles

struct DataClsp
    items::UnitRange{Int}
    periods::UnitRange{Int}
    numItems::Int
    numPer::Int
    cap::Int
    pc::Int
    sc
    hc
    pt
    st
    dem
end

struct Vertex
    vertex_id::Int64
    i::Int64            # item
    s::Int64            # inventory level
    t::Int64            # period
    inv_qtd::Int64      # inventory quantity
end

struct Arc
    arc_id::Int64
    i::Int64                    # item
    origin::Vertex              # origin vertex
    destination::Vertex         # destination vertex
    prod_qtd::Int64             # production amount of arc
    res_consumption::Int64      # resource consumption of arc
    cost::Float64               # cost of arc
end

struct InputGraph
    i::Int64            # item
    V::Array{Vertex}    # set of vertices
    A::Array{Arc}       # set of arcs
end

export DataClsp, read_clsp_data, compute_bigM_coeffs

function read_clsp_data(instanceFile)
   instance = readdlm(instanceFile)

   n = instance[1, 1]  # Get number of items
   items = 1:n
   m = instance[1, 2]  # Get number of periods
   periods = 1:m
   p = instance[2, 1]  # Get unitary production cost
   c = instance[3, 1]  # Get capacity

   a = Array{Float64}(undef, n)
   h = Array{Float64}(undef, n)
   b = Array{Float64}(undef, n)
   f = Array{Float64}(undef, n)
   for i=1:n
       a[i] = instance[3+i, 1]     # Get unitary resource consumptions
       h[i] = instance[3+i, 2]     # Get unitary inventory costs
       b[i] = instance[3+i, 3]     # Get setup resource consumptions
       f[i] = instance[3+i, 4]     # Get setup costs
   end

   d = Array{Int64}(undef, n, m)
   for i=1:n
       for t=1:m
           d[i, t] = instance[3+n+t, i]    # Get demands
       end
   end

   instance = DataClsp(items, periods, n, m, c, p, f, h, a, b, d)

   return instance

end

items(instance::DataClsp) = instance.items
periods(instance::DataClsp) = instance.periods
capacity(instance::DataClsp) = instance.cap
prod_cost(instance::DataClsp) = instance.pc
setup_cost(instance::DataClsp, i::Int) = instance.sc[i]
inv_cost(instance::DataClsp, i::Int) = instance.hc[i]
prod_time(instance::DataClsp, i::Int) = instance.pt[i]
setup_time(instance::DataClsp, i::Int) = instance.st[i]
demand(instance::DataClsp, i::Int, t::Int) = instance.dem[i, t]

function compute_bigM_coeffs(inst::DataClsp)
   bigM_coeffs = Array{Float64}(undef, inst.numItems, inst.numPer)

   for item = 1:inst.numItems
       for period = 1:inst.numPer
           bigM_coeffs[item, period] = (inst.cap - inst.st[item]) / inst.pt[item]

           temp = 0.0
           for t = period:inst.numPer
               temp += inst.dem[item, t]
           end

           if temp <= bigM_coeffs[item, period]
               bigM_coeffs[item, period] = temp
           end
       end
   end

   return bigM_coeffs
end

# Compute the sum of demands of item i from period t to period tt, with tt >= t
function sum_demands(inst::DataClsp, i::Int,  t::Int, tt::Int)
    sum = 0
    for period in t:tt
        sum += inst.dem[i, period]
    end
    return sum
end

# Build set S for item i
function build_WagnerWhitin_S(inst::DataClsp, i::Int)
    S = Array{Int}(undef, inst.numPer, inst.numPer + 1)

    fill!(S, -1)
    for period in 1:inst.numPer + 1
        S[1, period] = 0
    end

    for period in 2:inst.numPer
        max_inv_level = inst.numPer - period + 2
        for inv_level in 2:max_inv_level
            S[inv_level, period] = sum_demands(data, i, period, period + inv_level - 2)
        end
    end

    return S
end

# Build set of vertices for item i
function build_WagnerWhitin_V(inst::DataClsp, i::Int)
    S = build_WagnerWhitin_S(inst, i)
    V = Vertex[]
    id = 1
    for period in 1:inst.numPer + 1
        for inv_level in 1:inst.numPer
            if S[inv_level, period] != -1
                v = Vertex(id, i, inv_level, period, S[inv_level, period])
                push!(V, v)
                id += 1
            end
        end
    end

    return V
end

function build_S(inst::DataClsp, i::Int)
    S = Array{Int}(undef, sum_demands(inst, i, 1, inst.numPer), inst.numPer + 1)

    fill!(S, -1)
    for period in 1:inst.numPer + 1
        S[1, period] = 0
    end

    for period in 2:inst.numPer
        max_inv_level = min(floor((inst.cap - inst.st[i]) / inst.pt[i]), sum_demands(inst, i, 1, inst.numPer))
        for inv_level in 2:Int(max_inv_level)
            S[inv_level, period] = inv_level - 1
        end
    end

    return S
end

# Build set of vertices for item i
function build_V(inst::DataClsp, i::Int)
    S = build_S(inst, i)
    V = Vertex[]
    id = 1
    for period in 1:inst.numPer + 1
        for inv_level in 1:inst.numPer
            if S[inv_level, period] != -1
                v = Vertex(id, i, inv_level, period, S[inv_level, period])
                push!(V, v)
                id += 1
            end
        end
    end

    return V
end

function build_VV(inst::DataClsp, G::InputGraph)
    V = Int[]
    for v in G.V
        if v.t >= 2 && v.t <= inst.numPer
            push!(V, v.vertex_id)
        end
    end

    return V
end

# Build set of arcs for item i
function build_A(inst::DataClsp, i::Int, V::Array{Vertex})
    A = Arc[]

    id = 1
    for u in V
        for v in V
            if v.t == u.t + 1
                if inst.dem[i, u.t] - u.inv_qtd + v.inv_qtd >= 0
                    prod = inst.dem[i, u.t] - u.inv_qtd + v.inv_qtd
                    if prod == 0
                        cons = 0
                        cost = inst.hc[i] * v.inv_qtd
                    else
                        cons = (inst.pt[i] * prod) + inst.st[i]
                        cost = (inst.pc * prod) + (inst.hc[i] * v.inv_qtd) + inst.sc[i]
                    end
                    
                    a = Arc(id, i, u, v, prod, cons, cost)
                    push!(A, a)
                    id += 1
                end
            end
        end
    end

    return A
end

# Build the input graph for item i
function build_G(inst::DataClsp, i::Int)
    V = build_V(inst, i)
    A = build_A(inst, i, V)
    G = InputGraph(i, V, A)
    println("i = $i  G = ", G)

    return G
end

# return entering arcs in vertex v
function entering_arcs(G::InputGraph, v::Vertex)
    E = Int[]

    for a in G.A
        if a.destination == v
            push!(E, a.arc_id)
        end
    end

    return E
end

# return arcs leaving from vertex v
function leaving_arcs(G::InputGraph, v::Vertex)
    E = Int[]

    for a in G.A
        if a.origin == v
            push!(E, a.arc_id)
        end
    end

    return E
end


function show(io::IO, inst::DataClsp)
   # Print instance data
   println(io, "Instance data:")
   println(io, "Number of items: ", inst.numItems)
   println(io, "Number of periods: ", inst.numPer)
   println(io, "Unitary production cost: ", inst.pc)
   println(io, "Capacity in each period: ", inst.cap)
   println(io, "Items data:")
   println(io, "\tSetup costs: ", inst.sc)
   println(io, "\tUnitary inventory holding costs: ", inst.hc)
   println(io, "\tUnitary production times: ", inst.pt)
   println(io, "\tSetup resource consumptions: ", inst.st)
   for item = 1:inst.numItems   
      println(io, "Item $item demands: ", inst.dem[item, :])
   end

end
