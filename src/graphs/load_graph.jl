using Random
using Graphs
using SimpleWeightedGraphs
using JSON

# -------------------------------------------------------
# Detect Gset format
# -------------------------------------------------------
function _is_gset_file(path::AbstractString)
    open(path, "r") do f
        first = split(strip(readline(f)))
        return length(first) == 2 &&
               all(x -> tryparse(Int, x) !== nothing, first)
    end
end

# -------------------------------------------------------

function create_and_save_graph_QiILS(N::Int, k::Int, seed::Int;
    base_path::String = "/Users/guillermo.preisser/Projects/QiILS/graphs")

    # 1) generate unweighted regular graph
    g = Graphs.random_regular_graph(N, k; seed=seed)
    g_weighted = SimpleWeightedGraph(g)

    # 2) generate weights (same order as edges(g))
    Random.seed!(seed)
    nE = length(edges(g))
    rlist = rand(nE)

    # 3) assign weights in same order
    i = 1
    for e in edges(g)
        u, v = src(e), dst(e)
        w = rlist[i]
        Graphs.add_edge!(g_weighted, u, v, w)
        i += 1
    end

    # 4) Save to QiILS/graphs/N/k/
    dir_path = joinpath(base_path, string(N), string(k))
    mkpath(dir_path)

    filename = "graph_N$(N)_k$(k)_seed$(seed)_seedb$(seed).txt"
    full_path = joinpath(dir_path, filename)

    open(full_path, "w") do io
        for e in edges(g_weighted)
            u = src(e) - 1   # zero-based indexing to match Python
            v = dst(e) - 1
            w = g_weighted.weights[src(e), dst(e)]
            println(io, "$(w),$(u),$(v)")
        end
    end

    return g_weighted, full_path
end
# -------------------------------------------------------
# Custom CSV-like "<w,u,v>" graph format loader
# -------------------------------------------------------
function _load_custom_graph(path::AbstractString; weighted::Bool=true)
    lines = readlines(path)
    edges_parsed = [split(strip(ln), ",") for ln in lines]

    # infer number of nodes (0-based → shift to 1-based)
    nodes = Int[]
    for row in edges_parsed
        _, u, v = row
        push!(nodes, parse(Int, u) + 1)
        push!(nodes, parse(Int, v) + 1)
    end
    N = maximum(nodes)

    g = SimpleWeightedGraph(N)

    for row in edges_parsed
        w_str, u_str, v_str = row
        u = parse(Int, u_str) + 1
        v = parse(Int, v_str) + 1
        w = weighted ? parse(Float64, w_str) : 1.0
        add_edge!(g, u, v, w)
    end

    return g
end

# -------------------------------------------------------
# Gset loader (supports both 2- and 3-column formats)
# -------------------------------------------------------
function _load_gset(path::AbstractString)
    open(path, "r") do f
        first = split(strip(readline(f)))
        N = parse(Int, first[1])

        g = SimpleWeightedGraph(N)

        for line in eachline(f)
            parts = split(strip(line))

            if length(parts) == 2
                u = parse(Int, parts[1])
                v = parse(Int, parts[2])
                w = 1.0
            elseif length(parts) == 3
                u = parse(Int, parts[1])
                v = parse(Int, parts[2])
                w = parse(Float64, parts[3])
            else
                error("Invalid Gset line: $line")
            end

            add_edge!(g, u, v, w)
        end

        return g
    end
end

# -------------------------------------------------------
# Unified graph loader
# -------------------------------------------------------
function load_graph(; gset=nothing,
                     path=nothing,
                     N=nothing,
                     k=nothing,
                     weighted::Bool=true,
                     seed::Int=1)

    if gset !== nothing
        fname = "G$(gset)"
        if !isfile(fname)
            error("Gset file '$fname' not found.")
        end
        return _load_gset(fname)
    end

    if path !== nothing
        return _is_gset_file(path) ? _load_gset(path) :
                                     _load_custom_graph(path; weighted=weighted)
    end

    if N !== nothing && k !== nothing
        return _generate_random_regular(N, k; weighted=weighted, seed=seed)
    end

    error("load_graph requires gset=Int, path=String or (N,k).")
end