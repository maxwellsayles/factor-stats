set terminal eps enhanced

set xlabel 'Bits in Composite'

set key left

set xrange [50:62]
set ylabel 'Average Time (Microseconds)'
set output 'sspar-random-vs-sequential.eps'
plot 'sspar-factor/sspar-random-ideals-optimized-timings.dat' using 1:7 with lines t 'Random Prime Ideals', \
     'sspar-factor/sspar-optimized-timings.dat' using 1:7 with lines t 'Sequential Prime Ideals'
set xrange [*:*]

# Whiskers
set ylabel 'Time (Microseconds)'

set output 'factor-flint-whisker.eps'
plot 'flint-factor/flint-timings.dat' using 1:3:2:6:5 with candlesticks t 'Flint'

set output 'factor-pari-whisker.eps'
plot 'pari-factor/pari-timings.dat' using 1:3:2:6:5 with candlesticks t 'Pari'

set output 'factor-squfof-whisker.eps'
plot 'squfof-parigp/squfof-timings.dat' using 1:3:2:6:5 with candlesticks t 'SQUFOF'

set output 'factor-spar-whisker.eps'
plot 'spar-factor/spar-vanilla-timings.dat' using 1:3:2:6:5 with candlesticks t 'Vanilla SPAR'

set output 'factor-sspar-whisker.eps'
plot 'sspar-factor/sspar-optimized-timings.dat' using 1:3:2:6:5 with candlesticks t 'SuperSPAR'

# avg timings
set ylabel 'Average Time (Microseconds)'

set xrange [*:*]
set yrange [0:25000]
set output 'factor-average.eps'
plot 'flint-factor/flint-timings.dat' using 1:7 with lines t 'Flint', \
     'pari-factor/pari-timings.dat' using 1:7 with lines t 'Pari', \
     'squfof-parigp/squfof-timings.dat' using 1:7 with lines t 'SQUFOF', \
     'spar-factor/spar-vanilla-timings.dat' using 1:7 with lines t 'Vanilla SPAR', \
     'sspar-factor/sspar-optimized-timings.dat' using 1:7 with lines t 'SuperSPAR'
set yrange [*:*]

set xrange [46:58]
set output 'factor-average-zoom-left.eps'
plot 'flint-factor/flint-timings.dat' using 1:7 with lines t 'Flint', \
     'pari-factor/pari-timings.dat' using 1:7 with lines t 'Pari', \
     'squfof-parigp/squfof-timings.dat' using 1:7 with lines t 'SQUFOF', \
     'sspar-factor/sspar-optimized-timings.dat' using 1:7 with lines t 'SuperSPAR' lc 5

set xrange [58:66]
set output 'factor-average-zoom-right.eps'
plot 'flint-factor/flint-timings.dat' using 1:7 with lines t 'Flint', \
     'pari-factor/pari-timings.dat' using 1:7 with lines t 'Pari', \
     'squfof-parigp/squfof-timings.dat' using 1:7 with lines t 'SQUFOF', \
     'sspar-factor/sspar-optimized-timings.dat' using 1:7 with lines t 'SuperSPAR' lc 5

# median timings
set ylabel 'Median Time (Microseconds)'

set xrange [*:*]
set yrange [0:25000]
set output 'factor-median.eps'
plot 'flint-factor/flint-timings.dat' using 1:4 with lines t 'Flint', \
     'pari-factor/pari-timings.dat' using 1:4 with lines t 'Pari', \
     'squfof-parigp/squfof-timings.dat' using 1:4 with lines t 'SQUFOF', \
     'spar-factor/spar-vanilla-timings.dat' using 1:4 with lines t 'Vanilla SPAR', \
     'sspar-factor/sspar-optimized-timings.dat' using 1:4 with lines t 'SuperSPAR'
set yrange [*:*]

set xrange [44:58]
set output 'factor-median-zoom-left.eps'
plot 'flint-factor/flint-timings.dat' using 1:4 with lines t 'Flint', \
     'pari-factor/pari-timings.dat' using 1:4 with lines t 'Pari', \
     'squfof-parigp/squfof-timings.dat' using 1:4 with lines t 'SQUFOF', \
     'sspar-factor/sspar-optimized-timings.dat' using 1:4 with lines t 'SuperSPAR' lc 5

set xrange [58:70]
set output 'factor-median-zoom-right.eps'
plot 'flint-factor/flint-timings.dat' using 1:4 with lines t 'Flint', \
     'pari-factor/pari-timings.dat' using 1:4 with lines t 'Pari', \
     'squfof-parigp/squfof-timings.dat' using 1:4 with lines t 'SQUFOF', \
     'sspar-factor/sspar-optimized-timings.dat' using 1:4 with lines t 'SuperSPAR' lc 5