module model

using Parameters

export miRNA_model, rnaParamValues

function miRNA_model(du, u, p, t)
    mir339, mir486, mir6803, mir128, mir92, mir885, mir146a, mir125, mir504, mir146b = ones(10)
    NFKBT, BEC1T, CARD9T, GPR35T, INAVAT, GPR18T, APEHT, GPR65T, PARK7T = ones(10)

    # kanf=1.5 kanf1=1, kanf2=1, kabc=0.1, kabc1=2, kabc2=0.1
    #stress faktor
    kanf=0.15
    #kanf=1.5
    kanf1=0.1
    #kanf1=1.0
    kanf2=0.1
    #kanf2=1.0
    #
    kinf=1
    kinf1=1
    kinf2=15
    Jnfk=0.1
    #stress faktor
    kabc=0.01
    #kabc=0.1
    kabc1=0.2
    #kabc1=2.0
    kabc2=0.01
    #kabc2=0.1
    #
    kibc=0.5
    kibc1=1
    kibc2=15
    Jbec=0.1
    kacd9=1
    kicd9=0.1
    kagr35=1
    kigr35=0.1
    kaiaa=1
    kiiaa=0.1
    kaaph=1
    kiaph=0.1
    kagr18=1
    kigr18=0.1
    kagr65=1
    kigr65=0.1
    kapk7=1
    kipk7=0.1

    GPR1 = u[3]

    GPR2 = u[4] + u[5] + u[6]

    GPR3 = u[7] + u[8] + u[9]

    du[1] = (kanf + kanf1*GPR1 + kanf2*GPR2)*(NFKBT-u[1])/(Jnfk + NFKBT-u[1]) -
                (kinf + kinf1*GPR3 + kinf2*u[2])*u[1]/(Jnfk + u[1])

    du[2] = (kabc + kabc1*GPR2 + kabc2*GPR3)*(BEC1T-u[2])/(Jbec + BEC1T-u[2]) -
                (kibc + kibc1*GPR1 + kibc2*u[1])*u[2]/(Jbec + u[2])

    du[3] = kacd9*(CARD9T - u[3]) -
                (kicd9 + p.card9_mir6803*mir6803 + p.card9_mir486*mir486)*u[3]

    du[4] = kagr35*(GPR35T - u[4]) -
                (kigr35 + p.gpr35_mir6803*mir6803 + p.gpr35_mir128*mir128 +
                    p.gpr35_mir92*mir92 + p.gpr35_mir125*mir125)*u[4]

    du[5] = kaiaa*(INAVAT - u[5]) -
                (kiiaa + p.inava_mir504*mir504 + p.inava_mir146b*mir146b +
                p.inava_mir125*mir125)*u[5]

    du[6] = kaaph*(APEHT - u[6]) - (kiaph + p.apeh_mir504*mir504)*u[6]

    du[7] = kagr18*(GPR18T - u[7]) -
                (kigr18 + p.gpr18_mir885*mir885 + p.gpr18_mir146a*mir146a + p.gpr18_mir92*mir92 +
                p.gpr18_mir125*mir125 + p.gpr18_mir504*mir504 + p.gpr18_mir486*mir486)*u[7]

    du[8] = kagr65*(GPR65T - u[8]) -
                (kagr65 + kigr65 + p.gpr65_mir128*mir128)*u[8]

    du[9] = kapk7*(PARK7T - u[9]) -
                (kipk7 + p.park7_mir486*mir486 + p.park7_mir339*mir339)*u[9]
end
# stress esetén:
# kanf=1.5 kanf1=1, kanf2=1, kabc=0.1, kabc1=2, kabc2=0.1
end
