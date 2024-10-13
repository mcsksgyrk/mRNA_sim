module parsers

using XLSX
using DataFrames
using Parameters

export rnaParamValues, parseRNA, setParamValues

Base.@kwdef mutable struct rnaParamValues
    apeh_mir504::Float64 = 1.0
    card9_mir486::Float64 = 1.0
    card9_mir6803::Float64 = 0.0
    gpr18_mir125::Float64 = 1.0
    gpr18_mir146a::Float64 = 1.0
    gpr18_mir486::Float64 = 1.0
    gpr18_mir885::Float64 = 0.0
    gpr18_mir504::Float64 = 0.0
    gpr18_mir92::Float64 = 1.0
    gpr35_mir128::Float64 = 1.0
    gpr35_mir6803::Float64 = 1.0
    gpr35_mir92::Float64 = 1.0
    gpr35_mir125::Float64 = 1.0
    gpr65_mir128::Float64 = 1.0
    inava_mir125::Float64 = 0.0
    inava_mir146b::Float64 = 0.0
    inava_mir504::Float64 = 0.0
    park7_mir339::Float64 = 1.0
    park7_mir486::Float64 = 1.0
end
Base.copy(s::rnaParamValues) = rnaParamValues(s.apeh_mir504,s.card9_mir486, s.card9_mir6803,
                                              s.gpr18_mir125,s.gpr18_mir146a,s.gpr18_mir486,
                                              s.gpr18_mir885,s.gpr18_mir504,s.gpr18_mir92,
                                              s.gpr35_mir128,s.gpr35_mir6803,s.gpr35_mir92,
                                              s.gpr35_mir125,s.gpr65_mir128,s.inava_mir125,
                                              s.inava_mir146b,s.inava_mir504,s.park7_mir339,
                                              s.park7_mir486)
function parseRNA(row)
    mirna_id = "_mir"*split(row.miRNA,"-")[3]
    target_gene = row.Target_gene_name
    param_name = lowercase(target_gene*mirna_id)
    edge_value = abs(Float64(row.edge_weight))
    return param_name, edge_value
end

function setParamValues(paramSet,param_name,edge_value)
    try
        setfield!(paramSet,Symbol(param_name), edge_value)
    catch e
        sliced_name = param_name[1:end-1]
        println("no $param_name trying $sliced_name")
        setfield!(paramSet,Symbol(sliced_name), edge_value)
    end
end

end
