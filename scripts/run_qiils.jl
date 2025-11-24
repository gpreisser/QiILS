using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))   # activates QiILS root

using QiILS
using JSON
using Printf
using Graphs          # <-- ADD THIS
# ---------------------------------------------------------
# User parameters
# ---------------------------------------------------------
N        = 10
k        = 3
seed     = 2
weighted = true            # set false for unweighted graphs

λ_sweep = 1.0
attempts = 10
sweeps_per_attempt = 80
percentage = 0.1
angle_conv = 0.00000001

println("====================================================")
println("                QiILS Optimizer Runner              ")
println("====================================================")

# ---------------------------------------------------------
# Load graph
# ---------------------------------------------------------
println("▶ Loading graph N=$N, k=$k, seed=$seed, weighted=$weighted …")
wg, graphfile = create_and_save_graph_QiILS(N, k, seed)
println("Graph saved at $graphfile")
# ---------------------------------------------------------
# Load optimal cut (if available)
# ---------------------------------------------------------
solution_path = solution_file_path(N, k, seed; weighted=weighted)

optimal_cut = load_optimal_cut(solution_path)
if optimal_cut === nothing
    println("⚠ No optimal solution found. Will skip approximation ratio.")
else
    println("✔ Optimal MaxCut value = $optimal_cut")
end

# ---------------------------------------------------------
# Prepare gvec
# ---------------------------------------------------------
gvec = zeros(Float64, nv(wg))     # currently unused, safe placeholder

# ---------------------------------------------------------
# Run solver
# ---------------------------------------------------------
println("\n▶ Running QiILS solver…")
best_history, best_angles = qiils_solve(
    wg,
    λ_sweep,
    gvec;
    attempts=attempts,
    sweeps_per_attempt=sweeps_per_attempt,
    percentage=percentage,
    seed=seed,
    angle_conv=angle_conv,
)

best_cut = best_history[end]

println("\n====================================================")
println("                 QiILS Solver Results               ")
println("====================================================")
println("  ✓ Best MaxCut found = $best_cut")

if optimal_cut !== nothing
    ratio = best_cut / optimal_cut
    println("  ✓ Approximation ratio = $(round(ratio, digits=5))")
end

println("====================================================\n")

# ---------------------------------------------------------
# Save results
# ---------------------------------------------------------
weight_tag = weighted ? "weighted" : "unweighted"
save_dir = "results/N$(N)_k$(k)_seed$(seed)_$(weight_tag)/"
mkpath(save_dir)

save_data = Dict(
    "N" => N,
    "k" => k,
    "seed" => seed,
    "weighted" => weighted,
    "λ_sweep" => λ_sweep,
    "attempts" => attempts,
    "sweeps_per_attempt" => sweeps_per_attempt,
    "percentage" => percentage,
    "angle_conv" => angle_conv,
    "best_cut_history" => best_history,
    "best_cut" => best_cut,
    "best_angles" => best_angles,
    "optimal_cut" => optimal_cut === nothing ? "none" : optimal_cut,
)

save_path = save_dir * "qiils_results.json"
open(save_path, "w") do io
    JSON.print(io, save_data)
end

println("✔ Results saved to: $save_path\n")
println("Done.")