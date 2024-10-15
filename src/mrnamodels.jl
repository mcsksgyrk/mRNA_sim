include("src/submodules/model.jl")
include("src/submodules/parsers.jl")

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
data_df = DataFrame(XLSX.readtable("./source/patient_level_mirna_20240514.xlsx",
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
# basal state
basal_results = DataFrame([[],[],[],[],[]], ["patient", "NFKB_min", "NFKB_max", "BEC1_min", "BEC1_max"])
for (key, value) in param_sets
    p = value
    prob = ODEProblem(miRNA_model, u0, tspan, p)
    sol = solve(prob, Rosenbrock23())
    push!(basal_results, make_df_row(sol,key))
end
basal_results
############
stress_results = DataFrame([[],[],[],[],[]], ["patient", "NFKB_min", "NFKB_max", "BEC1_min", "BEC1_max"])
for (key, value) in param_sets
    p = value
    prob = ODEProblem(miRNA_model_stress, u0, tspan, p)
    sol = solve(prob, Rosenbrock23())
    push!(stress_results, make_df_row(sol,key))
end
maximum.(eachcol(stress_results))
minimum.(eachcol(stress_results))
maximum.(eachcol(basal_results))
#plot(sol.t,sol[1,:], ylim=(0, 1), label="NFKB")
#plot!(sol.t,sol[2,:], ylim=(0, 1), label="BEC1")
#sol.t
#stress_results[stress_results.patient="10252"]
## Extract the row labels (keys)
#rows = collect(keys(param_sets))  # ["Row1", "Row2", "Row3"]
#rows
## Dynamically extract all 19 fields from the struct into a matrix
#num_fields = fieldcount(rnaParamValues)  # Number of fields in the struct (19)
#field_names = fieldnames(rnaParamValues)
#data_matrix = [getfield(v, field) for v in values(param_sets), field in field_names]
## Plot heatmap
#heatmap(data_matrix,
#        xlabel="Fields",
#        ylabel="Rows",
#        xticks=(1:length(field_names), field_names),  # Dynamically label fields
#        yticks=(1:length(rows), rows),
#        color=:viridis,
#        title="Heatmap from Custom Structs with Unique Fields",
#        ticksrotation=90,
#        xrotation=90.0)
#
#fields_to_plot = [:card9_mir486, :gpr18_mir486]
#field_indices = findall(field -> field in fields_to_plot, field_names)
#heatmap(data_matrix,
#        xlabel="Fields",
#        ylabel="Patients",
#        xticks=(1:length(fields_to_plot), fields_to_plot),  # Label only the specified fields
#        yticks=(1:length(rows), rows),
#        color=:viridis,
#        title="Heatmap for Selected Fields",
#        xrotation=90)  # Rotate x-axis labels by 90 degrees
#using Clustering
#data_matrix
## Convert to a matrix
#data_matrix = hcat(data_matrix...)  # Transpose for correct orientation
#d = data_matrix+data_matrix'
## Perform hierarchical clustering
#dendrogram = hclust(data_matrix)  # Use transpose for correct dimensions
#
## Plot the dendrogram
#dendrogram_plot = dendrogram_plot(dendrogram, labels=keys(data_dict), title="Hierarchical Clustering Dendrogram")
#display(dendrogram_plot)
