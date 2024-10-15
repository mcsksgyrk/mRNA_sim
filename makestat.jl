using CSV
using DataFrames

basal = DataFrame(CSV.File("basal_results.csv"))
stress = DataFrame(CSV.File("stress_results.csv"))

basal.inflamated = basal.BEC1_max.<basal.NFKB_max
stress.inflamated = stress.BEC1_max.<stress.NFKB_max

CSV.write("stress_results_inf.csv",stress)
CSV.write("basal_results_inf.csv",basal)
