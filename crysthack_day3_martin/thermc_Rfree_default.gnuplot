set grid x y
set xrange [*:*] reverse
set xlabel "resolution (A)"
set ylabel "Rfree"
set key left center
set terminal png fontscale 2.0 size 900,700
set output 'thermc_Rfree_default.png'
plot './thermc_1-8A.csv' using 3:8 lw 2 ps 0.8 with linespoints title '1.8 A', './thermc_1-7A.csv' using 3:8 lw 2 ps 0.8 with linespoints title '1.7 A', './thermc_1-6A.csv' using 3:8 lw 2 ps 0.8 with linespoints title '1.6 A'
