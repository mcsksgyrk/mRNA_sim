using CSV
using DataFrames
using XLSX

basal = DataFrame(CSV.File("basal_results.csv"))
stress = DataFrame(CSV.File("stress_results.csv"))
patients = DataFrame(CSV.File("patient_configs.csv"))

basal.inflamated = basal.BEC1_max.<basal.NFKB_max
stress.inflamated = stress.BEC1_max.<stress.NFKB_max

inflamed_stress = filter(:inflamated => x -> x == true, stress)

inflamed_nop = filter(:inflamated => x -> x == false, stress)
nop_patients = inflamed_nop.patient
filtered_df2 = subset(patients, :ID => ByRow(x -> x in nop_patients))
CSV.write("bec1_active_configs.csv", filtered_df2)
#inava1 1 azok a sorok kellenek
inava1_ek = filter(:inava_mir125 => x -> x == 1.0, filtered_df2)
CSV.write("inava1_active_bec1_active.csv", inava1_ek)

kabumm_patients = inflamed_stress.patient
filtered_df = subset(patients, :ID => ByRow(x -> x in kabumm_patients))

#CSV.write("nfkb_active_configs.csv", filtered_df)
#CSV.write("stress_results_inf.csv",stress)
#CSV.write("basal_results_inf.csv",basal)
