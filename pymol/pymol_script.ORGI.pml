hide everything
viewport 800, 800
bg_color white
set ray_shadows,0
as sticks, all
set_bond stick_radius, 0.4, all

# for color by chromosome
spectrum b, rainbow, all, 1, 23

#################
# define groups #
#################

select genome, (b=0 or b=1 or b=2)
deselect

select nonOR, (b=0)
deselect

select GI, (b=1)
deselect

select OR, (b=2)
deselect

select mor28, (b=3)
deselect

select chr1_pat_nonOR,(b=0 and chain "chr1(pat)")
deselect
select chr1_mat_nonOR,(b=0 and chain "chr1(mat)")
deselect

select chr1_pat_GI,(b=1 and chain "chr1(pat)")
deselect
select chr1_mat_GI,(b=1 and chain "chr1(mat)")
deselect

select chr1_pat_OR,(b=2 and chain "chr1(pat)")
deselect
select chr1_mat_OR,(b=2 and chain "chr1(mat)")
deselect

select chr2_pat_nonOR,(b=0 and chain "chr2(pat)")
deselect
select chr2_mat_nonOR,(b=0 and chain "chr2(mat)")
deselect

select chr2_pat_GI,(b=1 and chain "chr2(pat)")
deselect
select chr2_mat_GI,(b=1 and chain "chr2(mat)")
deselect

select chr2_pat_OR,(b=2 and chain "chr2(pat)")
deselect
select chr2_mat_OR,(b=2 and chain "chr2(mat)")
deselect

select chr3_pat_nonOR,(b=0 and chain "chr3(pat)")
deselect
select chr3_mat_nonOR,(b=0 and chain "chr3(mat)")
deselect

select chr3_pat_GI,(b=1 and chain "chr3(pat)")
deselect
select chr3_mat_GI,(b=1 and chain "chr3(mat)")
deselect

select chr3_pat_OR,(b=2 and chain "chr3(pat)")
deselect
select chr3_mat_OR,(b=2 and chain "chr3(mat)")
deselect

select chr4_pat_nonOR,(b=0 and chain "chr4(pat)")
deselect
select chr4_mat_nonOR,(b=0 and chain "chr4(mat)")
deselect

select chr4_pat_GI,(b=1 and chain "chr4(pat)")
deselect
select chr4_mat_GI,(b=1 and chain "chr4(mat)")
deselect

select chr4_pat_OR,(b=2 and chain "chr4(pat)")
deselect
select chr4_mat_OR,(b=2 and chain "chr4(mat)")
deselect

select chr5_pat_nonOR,(b=0 and chain "chr5(pat)")
deselect
select chr5_mat_nonOR,(b=0 and chain "chr5(mat)")
deselect

select chr5_pat_GI,(b=1 and chain "chr5(pat)")
deselect
select chr5_mat_GI,(b=1 and chain "chr5(mat)")
deselect

select chr5_pat_OR,(b=2 and chain "chr5(pat)")
deselect
select chr5_mat_OR,(b=2 and chain "chr5(mat)")
deselect

select chr6_pat_nonOR,(b=0 and chain "chr6(pat)")
deselect
select chr6_mat_nonOR,(b=0 and chain "chr6(mat)")
deselect

select chr6_pat_GI,(b=1 and chain "chr6(pat)")
deselect
select chr6_mat_GI,(b=1 and chain "chr6(mat)")
deselect

select chr6_pat_OR,(b=2 and chain "chr6(pat)")
deselect
select chr6_mat_OR,(b=2 and chain "chr6(mat)")
deselect

select chr7_mat_nonOR,(b=0 and chain "chr7(mat)")
deselect
select chr7_pat_nonOR,(b=0 and chain "chr7(pat)")
deselect

select chr7_mat_GI,(b=1 and chain "chr7(mat)")
deselect
select chr7_pat_GI,(b=1 and chain "chr7(pat)")
deselect

select chr7_mat_OR,(b=2 and chain "chr7(mat)")
deselect
select chr7_pat_OR,(b=2 and chain "chr7(pat)")
deselect





select chr8_mat_nonOR,(b=0 and chain "chr8(mat)")
deselect
select chr8_pat_nonOR,(b=0 and chain "chr8(pat)")
deselect

select chr8_mat_GI,(b=1 and chain "chr8(mat)")
deselect
select chr8_pat_GI,(b=1 and chain "chr8(pat)")
deselect

select chr8_mat_OR,(b=2 and chain "chr8(mat)")
deselect
select chr8_pat_OR,(b=2 and chain "chr8(pat)")
deselect




select chr9_pat_nonOR,(b=0 and chain "chr9(pat)")
deselect
select chr9_mat_nonOR,(b=0 and chain "chr9(mat)")
deselect

select chr9_pat_GI,(b=1 and chain "chr9(pat)")
deselect
select chr9_mat_GI,(b=1 and chain "chr9(mat)")
deselect

select chr9_pat_OR,(b=2 and chain "chr9(pat)")
deselect
select chr9_mat_OR,(b=2 and chain "chr9(mat)")
deselect





select chr10_mat_nonOR,(b=0 and chain "chr10(mat)")
deselect
select chr10_pat_nonOR,(b=0 and chain "chr10(pat)")
deselect

select chr10_mat_GI,(b=1 and chain "chr10(mat)")
deselect
select chr10_pat_GI,(b=1 and chain "chr10(pat)")
deselect

select chr10_mat_OR,(b=2 and chain "chr10(mat)")
deselect
select chr10_pat_OR,(b=2 and chain "chr10(pat)")
deselect



select chr11_pat_nonOR,(b=0 and chain "chr11(pat)")
deselect
select chr11_mat_nonOR,(b=0 and chain "chr11(mat)")
deselect

select chr11_pat_GI,(b=1 and chain "chr11(pat)")
deselect
select chr11_mat_GI,(b=1 and chain "chr11(mat)")
deselect

select chr11_pat_OR,(b=2 and chain "chr11(pat)")
deselect
select chr11_mat_OR,(b=2 and chain "chr11(mat)")
deselect



select chr12_pat_nonOR,(b=0 and chain "chr12(pat)")
deselect
select chr12_mat_nonOR,(b=0 and chain "chr12(mat)")
deselect

select chr12_pat_GI,(b=1 and chain "chr12(pat)")
deselect
select chr12_mat_GI,(b=1 and chain "chr12(mat)")
deselect

select chr12_pat_OR,(b=2 and chain "chr12(pat)")
deselect
select chr12_mat_OR,(b=2 and chain "chr12(mat)")
deselect



select chr13_pat_nonOR,(b=0 and chain "chr13(pat)")
deselect
select chr13_mat_nonOR,(b=0 and chain "chr13(mat)")
deselect

select chr13_pat_GI,(b=1 and chain "chr13(pat)")
deselect
select chr13_mat_GI,(b=1 and chain "chr13(mat)")
deselect

select chr13_pat_OR,(b=2 and chain "chr13(pat)")
deselect
select chr13_mat_OR,(b=2 and chain "chr13(mat)")
deselect



select chr14_mat_nonOR,(b=0 and chain "chr14(mat)")
deselect
select chr14_pat_nonOR,(b=0 and chain "chr14(pat)")
deselect

select chr14_mat_GI,(b=1 and chain "chr14(mat)")
deselect
select chr14_pat_GI,(b=1 and chain "chr14(pat)")
deselect

select chr14_mat_OR,(b=2 and chain "chr14(mat)")
deselect
select chr14_pat_OR,(b=2 and chain "chr14(pat)")
deselect

select chr14_mat_mor28,(b=3 and chain "chr14(mat)")
deselect
select chr14_pat_mor28,(b=3 and chain "chr14(pat)")
deselect



select chr15_mat_nonOR,(b=0 and chain "chr15(mat)")
deselect
select chr15_pat_nonOR,(b=0 and chain "chr15(pat)")
deselect

select chr15_mat_GI,(b=1 and chain "chr15(mat)")
deselect
select chr15_pat_GI,(b=1 and chain "chr15(pat)")
deselect

select chr15_mat_OR,(b=2 and chain "chr15(mat)")
deselect
select chr15_pat_OR,(b=2 and chain "chr15(pat)")
deselect



select chr16_pat_nonOR,(b=0 and chain "chr16(pat)")
deselect
select chr16_mat_nonOR,(b=0 and chain "chr16(mat)")
deselect

select chr16_pat_GI,(b=1 and chain "chr16(pat)")
deselect
select chr16_mat_GI,(b=1 and chain "chr16(mat)")
deselect

select chr16_pat_OR,(b=2 and chain "chr16(pat)")
deselect
select chr16_mat_OR,(b=2 and chain "chr16(mat)")
deselect



select chr17_mat_nonOR,(b=0 and chain "chr17(mat)")
deselect
select chr17_pat_nonOR,(b=0 and chain "chr17(pat)")
deselect

select chr17_mat_GI,(b=1 and chain "chr17(mat)")
deselect
select chr17_pat_GI,(b=1 and chain "chr17(pat)")
deselect

select chr17_mat_OR,(b=2 and chain "chr17(mat)")
deselect
select chr17_pat_OR,(b=2 and chain "chr17(pat)")
deselect



select chr18_mat_nonOR,(b=0 and chain "chr18(mat)")
deselect
select chr18_pat_nonOR,(b=0 and chain "chr18(pat)")
deselect

select chr18_mat_GI,(b=1 and chain "chr18(mat)")
deselect
select chr18_pat_GI,(b=1 and chain "chr18(pat)")
deselect

select chr18_mat_OR,(b=2 and chain "chr18(mat)")
deselect
select chr18_pat_OR,(b=2 and chain "chr18(pat)")
deselect



select chr19_pat_nonOR,(b=0 and chain "chr19(pat)")
deselect
select chr19_mat_nonOR,(b=0 and chain "chr19(mat)")
deselect

select chr19_pat_GI,(b=1 and chain "chr19(pat)")
deselect
select chr19_mat_GI,(b=1 and chain "chr19(mat)")
deselect

select chr19_pat_OR,(b=2 and chain "chr19(pat)")
deselect
select chr19_mat_OR,(b=2 and chain "chr19(mat)")
deselect



select chr20_pat_nonOR,(b=0 and chain "chr20(pat)")
deselect
select chr20_mat_nonOR,(b=0 and chain "chr20(mat)")
deselect

select chr20_pat_GI,(b=1 and chain "chr20(pat)")
deselect
select chr20_mat_GI,(b=1 and chain "chr20(mat)")
deselect

select chr20_pat_OR,(b=2 and chain "chr20(pat)")
deselect
select chr20_mat_OR,(b=2 and chain "chr20(mat)")
deselect



select chrX_mat_nonOR,(b=0 and chain "chrX(mat)")
deselect
select chrX_pat_nonOR,(b=0 and chain "chrX(pat)")
deselect

select chrX_mat_GI,(b=1 and chain "chrX(mat)")
deselect
select chrX_pat_GI,(b=1 and chain "chrX(pat)")
deselect

select chrX_mat_OR,(b=2 and chain "chrX(mat)")
deselect
select chrX_pat_OR,(b=2 and chain "chrX(pat)")
deselect


##################
# define colors #
#################

set_color chr1_color, [230, 25, 75]
set_color chr2_color, [60, 180, 75]
set_color chr3_color, [255, 225, 25]
set_color chr4_color, [67, 99, 216]
set_color chr5_color, [245, 130, 49]
set_color chr6_color, [145, 30, 180]
set_color chr7_color, [66, 212, 244]
set_color chr8_color, [240, 50, 230]
set_color chr9_color, [191, 239, 69]
set_color chr10_color, [250, 190, 190]
set_color chr11_color, [70, 153, 144]
set_color chr12_color, [230, 190, 255]
set_color chr13_color, [154, 99, 36]
set_color chr14_color, [255, 250, 200]
set_color chr15_color, [128, 0, 0]
set_color chr16_color, [170, 255, 195]
set_color chr17_color, [128, 128, 0]
set_color chr18_color, [255, 216, 177]
set_color chr19_color, [0, 0, 117]
set_color chr20_color, [169, 169, 169]
set_color chrX_color, [255, 255, 255]
set_color chrY_color, [0, 0, 0]
set_color nonOR_color, [210,210,210]

#########################################
# color things  for zoomed out version  #
#########################################

## NOTE: Setting -> Transparency -> Multi-Layer (Real Time OIT)
#
#hide everything, genome
#show sticks, genome
#set stick_transparency, 0.98, genome
#color nonOR_color, (genome)
#
## for cell 30 and cell 22
#
#show sticks, (chr9_pat_nonOR)
#show sticks, (chr9_pat_GI)
#show sticks, (chr9_pat_OR)
#set stick_transparency, 0.2, chr9_pat_nonOR
#set stick_transparency, 0.2, chr9_pat_GI
#set stick_transparency, 0.2, chr9_pat_OR
#color chr9_color, (chr9_pat_nonOR)
#color chr9_color, (chr9_pat_GI)
#color chr9_color, (chr9_pat_OR)
#show spheres, chr9_pat_OR
#set sphere_transparency, 0.0, chr9_pat_OR
#set sphere_scale, 0.5, chr9_pat_OR
#
#show sticks, (chr19_pat_nonOR)
#show sticks, (chr19_pat_GI)
#show sticks, (chr19_pat_OR)
#set stick_transparency, 0.4, chr19_pat_nonOR
#set stick_transparency, 0.4, chr19_pat_GI
#set stick_transparency, 0.4, chr19_pat_OR
#color chr19_color, (chr19_pat_nonOR)
#color chr19_color, (chr19_pat_GI)
#color chr19_color, (chr19_pat_OR)
#show spheres, chr19_pat_OR
#set sphere_transparency, 0.0, chr19_pat_OR
#set sphere_scale, 0.5, chr19_pat_OR
#
#show sticks, (chr2_mat_nonOR)
#show sticks, (chr2_mat_GI)
#show sticks, (chr2_mat_OR)
#set stick_transparency, 0.2, chr2_mat_nonOR
#set stick_transparency, 0.2, chr2_mat_GI
#set stick_transparency, 0.2, chr2_mat_OR
#color chr2_color, (chr2_mat_nonOR)
#color chr2_color, (chr2_mat_GI)
#color chr2_color, (chr2_mat_OR)
#show spheres, chr2_mat_OR
#set sphere_transparency, 0.0, chr2_mat_OR
#set sphere_scale, 0.5, chr2_mat_OR
#
## for cell 13
#
##show sticks, (chr9_pat_nonOR)
##show sticks, (chr9_pat_GI)
##show sticks, (chr9_pat_OR)
##set stick_transparency, 0.2, chr9_pat_nonOR
##set stick_transparency, 0.2, chr9_pat_GI
##set stick_transparency, 0.2, chr9_pat_OR
##color chr9_color, (chr9_pat_nonOR)
##color chr9_color, (chr9_pat_GI)
##color chr9_color, (chr9_pat_OR)
##show spheres, chr9_pat_OR
##set sphere_transparency, 0.0, chr9_pat_OR
##
##show sticks, (chr19_mat_nonOR)
##show sticks, (chr19_mat_GI)
##show sticks, (chr19_mat_OR)
##set stick_transparency, 0.4, chr19_mat_nonOR
##set stick_transparency, 0.4, chr19_mat_GI
##set stick_transparency, 0.4, chr19_mat_OR
##color chr19_color, (chr19_mat_nonOR)
##color chr19_color, (chr19_mat_GI)
##color chr19_color, (chr19_mat_OR)
##show spheres, chr19_mat_OR
##set sphere_transparency, 0.0, chr19_mat_OR
##
##show sticks, (chr2_pat_nonOR)
##show sticks, (chr2_pat_GI)
##show sticks, (chr2_pat_OR)
##set stick_transparency, 0.2, chr2_pat_nonOR
##set stick_transparency, 0.2, chr2_pat_GI
##set stick_transparency, 0.2, chr2_pat_OR
##color chr2_color, (chr2_pat_nonOR)
##color chr2_color, (chr2_pat_GI)
##color chr2_color, (chr2_pat_OR)
##show spheres, chr2_pat_OR
##set sphere_transparency, 0.0, chr2_pat_OR
#
##set sphere_scale, 0.5
#
##png z5OMP.13.chr9_chr19_chr2.whole_chr.png
##png z5OMP.22.chr9_chr19_chr2.whole_chr.png

####################
# ZOOMED IN VERSION #
#####################

#hide everything, genome
#
##set_bond stick_radius, 0.2, (chr19_pat_nonOR or chr19_pat_OR or chr19_pat_GI)
##set_bond stick_radius, 0.2, (chr2_mat_nonOR or chr2_mat_OR or chr2_mat_GI)
##set_bond stick_radius, 0.2, (chr9_pat_nonOR or chr9_pat_OR or chr9_pat_GI)
#
## for cell 30 and 22
#show spheres, (chr19_pat_OR or chr9_pat_OR or chr2_mat_OR)
#set sphere_scale, 0.2, (chr19_pat_OR or chr9_pat_OR or chr2_mat_OR)
#
## for cell 13
##show spheres, (chr19_mat_OR or chr9_pat_OR or chr2_pat_OR)
##set sphere_scale, 0.2, (chr19_mat_OR or chr9_pat_OR or chr2_pat_OR)
#
## all GIs
#show spheres, GI
#set sphere_scale, 0.7, GI
#set sphere_transparency, 0.00, GI
#
#color chr1_color, (chr1_pat_GI)
#color chr1_color, (chr1_mat_GI)
#color chr2_color, (chr2_pat_GI)
#color chr2_color, (chr2_mat_GI)
#color chr3_color, (chr3_pat_GI)
#color chr3_color, (chr3_mat_GI)
#color chr4_color, (chr4_pat_GI)
#color chr4_color, (chr4_mat_GI)
#color chr5_color, (chr5_pat_GI)
#color chr5_color, (chr5_mat_GI)
#color chr6_color, (chr6_pat_GI)
#color chr6_color, (chr6_mat_GI)
#color chr7_color, (chr7_pat_GI)
#color chr7_color, (chr7_mat_GI)
#color chr8_color, (chr8_pat_GI)
#color chr8_color, (chr8_mat_GI)
#color chr9_color, (chr9_pat_GI)
#color chr9_color, (chr9_mat_GI)
#color chr10_color, (chr10_pat_GI)
#color chr10_color, (chr10_mat_GI)
#color chr11_color, (chr11_pat_GI)
#color chr11_color, (chr11_mat_GI)
#color chr12_color, (chr12_pat_GI)
#color chr12_color, (chr12_mat_GI)
#color chr13_color, (chr13_pat_GI)
#color chr13_color, (chr13_mat_GI)
#color chr14_color, (chr14_pat_GI)
#color chr14_color, (chr14_mat_GI)
#color chr15_color, (chr15_pat_GI)
#color chr15_color, (chr15_mat_GI)
#color chr16_color, (chr16_pat_GI)
#color chr16_color, (chr16_mat_GI)
#color chr17_color, (chr17_pat_GI)
#color chr17_color, (chr17_mat_GI)
#color chr18_color, (chr18_pat_GI)
#color chr18_color, (chr18_mat_GI)
#color chr19_color, (chr19_pat_GI)
#color chr19_color, (chr19_mat_GI)
#color chr20_color, (chr20_pat_GI)
#color chr20_color, (chr20_mat_GI)
#color chrX_color, (chrX_pat_GI)
#color chrX_color, (chrX_mat_GI)
#
#set fog_start, 0.3
#set fog, 2
#
#png z5OMP.22.chr9_chr19_chr2.zoomed.png

################################
# color all chromoeoms version #
################################

hide everything, nonOR
hide everything, OR
show spheres, OR
set sphere_scale, 0.2, OR
hide everything, GI
show spheres, GI
set sphere_scale, 0.5, GI
hide everything, chr14_pat_mor28
show spheres, chr14_pat_mor28
set sphere_scale, 1.5, chr14_pat_mor28



color chr1_color, (chr1_pat_GI, chr1_pat_OR, chr1_pat_nonOR)
color chr1_color, (chr1_mat_GI, chr1_mat_OR, chr1_mat_nonOR)
color chr2_color, (chr2_pat_GI, chr2_pat_OR, chr2_pat_nonOR)
color chr2_color, (chr2_mat_GI, chr2_mat_OR, chr2_mat_nonOR)
color chr3_color, (chr3_pat_GI, chr3_pat_OR, chr3_pat_nonOR)
color chr3_color, (chr3_mat_GI, chr3_mat_OR, chr3_mat_nonOR)
color chr4_color, (chr4_pat_GI, chr4_pat_OR, chr4_pat_nonOR)
color chr4_color, (chr4_mat_GI, chr4_mat_OR, chr4_mat_nonOR)
color chr5_color, (chr5_pat_GI, chr5_pat_OR, chr5_pat_nonOR)
color chr5_color, (chr5_mat_GI, chr5_mat_OR, chr5_mat_nonOR)
color chr6_color, (chr6_pat_GI, chr6_pat_OR, chr6_pat_nonOR)
color chr6_color, (chr6_mat_GI, chr6_mat_OR, chr6_mat_nonOR)
color chr7_color, (chr7_pat_GI, chr7_pat_OR, chr7_pat_nonOR)
color chr7_color, (chr7_mat_GI, chr7_mat_OR, chr7_mat_nonOR)
color chr8_color, (chr8_pat_GI, chr8_pat_OR, chr8_pat_nonOR)
color chr8_color, (chr8_mat_GI, chr8_mat_OR, chr8_mat_nonOR)
color chr9_color, (chr9_pat_GI, chr9_pat_OR, chr9_pat_nonOR)
color chr9_color, (chr9_mat_GI, chr9_mat_OR, chr9_mat_nonOR)
color chr10_color, (chr10_pat_GI, chr10_pat_OR, chr10_pat_nonOR)
color chr10_color, (chr10_mat_GI, chr10_mat_OR, chr10_mat_nonOR)
color chr11_color, (chr11_pat_GI, chr11_pat_OR, chr11_pat_nonOR)
color chr11_color, (chr11_mat_GI, chr11_mat_OR, chr11_mat_nonOR)
color chr12_color, (chr12_pat_GI, chr12_pat_OR, chr12_pat_nonOR)
color chr12_color, (chr12_mat_GI, chr12_mat_OR, chr12_mat_nonOR)
color chr13_color, (chr13_pat_GI, chr13_pat_OR, chr13_pat_nonOR)
color chr13_color, (chr13_mat_GI, chr13_mat_OR, chr13_mat_nonOR)
color chr14_color, (chr14_pat_GI, chr14_pat_OR, chr14_pat_nonOR, chr14_pat_mor28)
color chr14_color, (chr14_mat_GI, chr14_mat_OR, chr14_mat_nonOR, chr14_mat_mor28)
color chr15_color, (chr15_pat_GI, chr15_pat_OR, chr15_pat_nonOR)
color chr15_color, (chr15_mat_GI, chr15_mat_OR, chr15_mat_nonOR)
color chr16_color, (chr16_pat_GI, chr16_pat_OR, chr16_pat_nonOR)
color chr16_color, (chr16_mat_GI, chr16_mat_OR, chr16_mat_nonOR)
color chr17_color, (chr17_pat_GI, chr17_pat_OR, chr17_pat_nonOR)
color chr17_color, (chr17_mat_GI, chr17_mat_OR, chr17_mat_nonOR)
color chr18_color, (chr18_pat_GI, chr18_pat_OR, chr18_pat_nonOR)
color chr18_color, (chr18_mat_GI, chr18_mat_OR, chr18_mat_nonOR)
color chr19_color, (chr19_pat_GI, chr19_pat_OR, chr19_pat_nonOR)
color chr19_color, (chr19_mat_GI, chr19_mat_OR, chr19_mat_nonOR)
color chr20_color, (chr20_pat_GI, chr20_pat_OR, chr20_pat_nonOR)
color chr20_color, (chr20_mat_GI, chr20_mat_OR, chr20_mat_nonOR)
color chrX_color, (chrX_pat_GI, chrX_pat_OR, chrX_pat_nonOR)
color chrX_color, (chrX_mat_GI, chrX_mat_OR, chrX_mat_nonOR)



hide everything, nonOR
hide everything, OR
show spheres, OR
set sphere_scale, 0.2, OR
set sphere_transparency, 0, OR
hide everything, GI
show spheres, GI
set sphere_scale, 0.5, GI
hide everything, chr14_pat_mor28
show spheres, chr14_pat_mor28
set sphere_scale, 1.5, chr14_pat_mor28


hide everything, nonOR
hide everything, OR
show spheres, OR
set sphere_scale, 0.2, OR
set sphere_transparency, 0.8, OR
hide everything, GI
show spheres, GI
set sphere_scale, 0.5, GI
hide everything, chr14_pat_mor28
show spheres, chr14_pat_mor28
set sphere_scale, 1.5, chr14_pat_mor28



hide surface,nonOR
show spheres,nonOR
set sphere_scale, 0.5,nonOR
set sphere_transparency, 0, nonOR
set sphere_transparency, 0, OR
set sphere_transparency, 0, GI


mset 1x60
movie_fade sphere_transparency, 1, 0, 60, 1, nonOR

mclear
mset 1 x360
frame 1
mview store
frame 60
movie_fade sphere_transparency, 1, 0, 60, 1, nonOR
mview store
turn y, 60
mview store
mplay

frame 180
zoom chr14_pat_OR
mview store
frame 300
zoom chr14_pat_OR
turn y, 120
mview store
mview interpolate
mplay


