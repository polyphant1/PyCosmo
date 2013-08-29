


norm = evs.void_survey_norm(vd,cosm)

pdf = []

for i in radius:
    pdf.append(evs.void_survey_pdf(i,norm,vd,10**5))
