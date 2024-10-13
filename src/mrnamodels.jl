include("submodules/model.jl")
include("submodules/parsers.jl")

using .model
using .parsers
using DifferentialEquations
using Plots
using XLSX
using DataFrames

function make_df_row(sol, patient)
    return [patient,
            minimum(sol[1,:]),
            maximum(sol[1,:]),
            minimum(sol[2,:]),
            maximum(sol[2,:])]
end

# Parse patient data
data_df = DataFrame(XLSX.readtable("../source/patient_level_mirna_20240514.xlsx",
                                   "patient_configurations"))
sort!(data_df, :Patient)
patients = sort(unique(values(data_df.Patient)))
rnaParams = rnaParamValues()

param_sets = Dict{String, rnaParamValues}()
for p in patients
    df = data_df[(data_df.Patient.==p), :]
    p_params = copy(rnaParams)
    for row in eachrow(df)
        param_name, edge_value = parseRNA(row)
        setParamValues(p_params, param_name, edge_value)
    end
    push!(param_sets,p=>p_params)
end
#########################################
u0 =[0.01601;
    0.03322;
    0.47619;
    0.19608;
    0.90909;
    0.47619;
    0.19608;
    0.32258;
    0.32258]
tspan = (0.0,100.0)
stress_free_results = DataFrame([[],[],[],[],[]], ["patient", "BEC1_min", "BEC1_max", "NFKB_min", "NFKB_max"])
for (key, value) in param_sets
    p = value
    prob = ODEProblem(miRNA_model, u0, tspan, p)
    sol = solve(prob, Rosenbrock23())
    push!(stress_free_results, make_df_row(sol,key))
end
stress_free_results
plot(sol.t,sol[1,:], ylim=(0, 1), label="BEC1")
plot!(sol.t,sol[2,:], ylim=(0, 1), label="NFKB")
sol.t

