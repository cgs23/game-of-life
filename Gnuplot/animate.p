set terminal gif animate enhanced delay 15 font "Times" 20 size 1280, 720
set output "gol.gif"
unset key
unset border
unset xtics
unset ytics
set size ratio -1
unset grid
gens = system("find . -name '*.txt' | wc  -l") 
gen(i) = sprintf("Generations: %5d\n", i)
initfile = "generation1.txt"
stats initfile nooutput
rows = cols = int(sqrt(STATS_records))
cells = STATS_records
array a[cells]
set xr[0:rows]
set yr[cols:0]
do for [g=1:gens]{
idx = 0
pos = 1
do for [i=0:rows-1]{
col = 1
do for [b=0:sqrt(blocks)-1]{
do for [j=1:sqrt(cells/blocks)] {
m = idx + b*int(cells/blocks) + j
a[pos] = int(system("cat generation".g.".txt | awk ''NR==". m ."''"))
if(a[pos]==0){
set object pos rectangle from col-1, i to col, i+1 lw 1 fc rgb "white"
} else {
set object pos rectangle from col-1, i to col, i+1 lw 1 fc rgb "black"
}
pos = pos + 1
col = col + 1
}
}
if((i+1)%int(sqrt(cells/blocks))==0){
idx = idx + int(sqrt(cells/blocks)) + int(cells/blocks)
} else{
idx = idx + int(sqrt(cells/blocks))
}
}
set label 1 center gen(g) font 'Verdana, 20' at graph 0.50, 1.015
plot 1/0
}
set out