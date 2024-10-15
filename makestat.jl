using CSV
using DataFrames
using XLSX

basal = DataFrame(CSV.File("basal_results.csv"))
stress = DataFrame(CSV.File("stress_results.csv"))
patients = DataFrame(CSV.File("patient_configs.csv"))

basal.inflamated = basal.BEC1_max.<basal.NFKB_max
stress.inflamated = stress.BEC1_max.<stress.NFKB_max

inflamed_stress = filter(:inflamated => x -> x == true, stress)

kabumm_patients = inflamed_stress.patient
filtered_df = subset(patients, :ID => ByRow(x -> x in kabumm_patients))

CSV.write("nfkb_active_configs.csv", filtered_df)

#CSV.write("stress_results_inf.csv",stress)
#CSV.write("basal_results_inf.csv",basal)
