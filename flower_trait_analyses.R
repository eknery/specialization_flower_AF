ftraits = read.table("0_data/flower_trait_matrix.csv", sep=",", h=T)


### choosing diagnostic traits
# select usefull traits
ftraits = ftraits[,c(1,3,4,8,10,12,14,16,17,18)]
str(ftraits)

# compound traits
major_stamen_height = ftraits$major_filet_height + ftraits$major_anther_height
minor_stamen_height = ftraits$minor_filet_height + ftraits$minor_anther_height
stamen_dimetrism = major_stamen_height - minor_stamen_height
herkogamy_major = ftraits$style_height - major_stamen_height
herkogamy_minor = ftraits$style_height - minor_stamen_height

 plot(stamen_dimetrism, ftraits$pore_long_section)
