set terminal eps enhanced

set xrange 'Bits in Composite'
set yrange 'Average Time (Microseconds)'
set key left

set xrange [50:62]
set output 'sspar-random-vs-sequential.eps'
plot 'sspar-random-ideals-optimized-timings.dat' using 1:7 with lines t 'Random Prime Ideals', \
     'sspar-optimized-timings.dat' using 1:7 with lines t 'Sequential Prime Ideals'
