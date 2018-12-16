##########################################################################################
# This data was published by Hannah E. Marx et. al. (2015), and need to be cited if used.
##########################################################################################

require(data.table)

traits <- "https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fddi.12401&file=ddi12401-sup-0012-AppendixS2.csv"

incidence_matrix <- "https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fddi.12401&file=ddi12401-sup-0013-AppendixS3.csv"

traits <- fread(traits)
incidence_matrix <- fread(incidence_matrix)
