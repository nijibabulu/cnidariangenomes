set terminal postscript color solid "Helvetica" 15
set output "figures/EDF3_alignment.ps"
set xtics rotate ( \
 "chr1" 1.0, \
 "chr2" 21651156.0, \
 "chr3" 43119771.0, \
 "chr4" 63945243.0, \
 "chr5" 82264127.0, \
 "chr6" 99965201.0, \
 "chr7" 116912359.0, \
 "chr8" 133644087.0, \
 "chr9" 149814405.0, \
 "chr10" 164999445.0, \
 "chr11" 179919353.0, \
 "chr12" 194716863.0, \
 "chr13" 208911330.0, \
 "chr14" 222708686.0, \
 "chr15" 235836516.0, \
 "" 247693854.0, \
 "" 250032310.0, \
 "" 250670031.0, \
 "" 251113484.0, \
 "" 251286479.0, \
 "" 251422026.0, \
 "" 251549798.0, \
 "" 251668039.0, \
 "" 251785739.0, \
 "" 251877413.0, \
 "" 251965710.0, \
 "" 252042083.0, \
 "" 252108937.0, \
 "" 252160337.0, \
 "" 252191773.0, \
 "" 252196800 \
)
set ytics ( \
 "chr1" 1.0, \
 "chr2" 21964434.0, \
 "chr3" 43249769.0, \
 "chr4" 64215161.0, \
 "chr5" 82687405.0, \
 "chr6" 100189483.0, \
 "chr7" 117168047.0, \
 "chr8" 134068912.0, \
 "chr9" 150845935.0, \
 "chr10" 166081305.0, \
 "chr11" 181024518.0, \
 "chr12" 195724775.0, \
 "chr13" 209699694.0, \
 "chr14" 222869273.0, \
 "chr15" 235994807.0, \
 "" 247797502.0, \
 "" 249694926.0, \
 "" 250307908.0, \
 "" 250751361.0, \
 "" 250963802.0, \
 "" 251137137.0, \
 "" 251260552.0, \
 "" 251367686.0, \
 "" 251442535.0, \
 "" 251500420.0, \
 "" 251550443.0, \
 "" 251593811.0, \
 "" 251630093.0, \
 "" 251664278.0, \
 "" 251695094.0, \
 "" 251722792.0, \
 "" 251749874.0, \
 "" 251776346.0, \
 "" 251801602.0, \
 "" 251826106.0, \
 "" 251849925.0, \
 "" 251872967.0, \
 "" 251894897.0, \
 "" 251916539.0, \
 "" 251937721.0, \
 "" 251957997.0, \
 "" 251978187.0, \
 "" 251997568.0, \
 "" 252013260.0, \
 "" 252028756.0, \
 "" 252043467.0, \
 "" 252056602.0, \
 "" 252068709.0, \
 "" 252079982.0, \
 "" 252090537.0, \
 "" 252100903.0, \
 "" 252110826.0, \
 "" 252120626.0, \
 "" 252129336.0, \
 "" 252134847.0, \
 "" 252138893.0, \
 "" 252142891.0, \
 "" 252146635.0, \
 "" 252148984.0, \
 "" 252151229.0, \
 "" 252153390.0, \
 "" 252155467.0, \
 "" 252157431.0, \
 "" 252158600 \
)
set size 1,1
set grid
unset key
set border 0
set tics scale 0
set xlabel "nv_contigs"
set ylabel "nv_dovetail"
set format "%.0f"
set mouse format "%.0f"
set mouse mouseformat "[%.0f, %.0f]"
if(GPVAL_VERSION < 5) { set mouse clipboardformat "[%.0f, %.0f]" }
set xrange [1:252196800]
set yrange [1:252158600]
set style line 1  lt 1 lw 2 pt 6 ps 0.25
set style line 2  lt 3 lw 2 pt 6 ps 0.25
set style line 3  lt 2 lw 2 pt 6 ps 0.25
plot \
 "data/assemblies/alignments/nv_dovetail-nv_rb.out.fplot" title "FWD" w lp ls 1, \
 "data/assemblies/alignments/nv_dovetail-nv_rb.out.rplot" title "REV" w lp ls 2
