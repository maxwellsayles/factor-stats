set terminal eps enhanced

set xlabel 'Bits in Composite'
set ylabel 'Average Time (Microseconds)'
set key left

set xrange [50:62]
set output 'sspar-random-vs-sequential.eps'
plot 'sspar-factor/sspar-random-ideals-optimized-timings.dat' using 1:7 with lines t 'Random Prime Ideals', \
     'sspar-factor/sspar-optimized-timings.dat' using 1:7 with lines t 'Sequential Prime Ideals'
set xrange [*:*]

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
set output 'factor-average.eps'
plot 'flint-factor/flint-timings.dat' using 1:7 with lines t 'Flint', \
     'pari-factor/pari-timings.dat' using 1:7 with lines t 'Pari', \
     'squfof-parigp/squfof-timings.dat' using 1:7 with lines t 'SQUFOF', \
     'spar-factor/spar-vanilla-timings.dat' using 1:7 with lines t 'Vanilla SPAR', \
     'sspar-factor/sspar-optimized-timings.dat' using 1:7 with lines t 'SuperSPAR'

set xrange [46:58]
set output 'factor-average-zoom-left.eps'
plot 'flint-factor/flint-timings.dat' using 1:7 with lines t 'Flint', \
     'pari-factor/pari-timings.dat' using 1:7 with lines t 'Pari', \
     'squfof-parigp/squfof-timings.dat' using 1:7 with lines t 'SQUFOF', \
     'spar-factor/spar-vanilla-timings.dat' using 1:7 with lines t 'Vanilla SPAR', \
     'sspar-factor/sspar-optimized-timings.dat' using 1:7 with lines t 'SuperSPAR'

set xrange [58:66]
set output 'factor-average-zoom-right.eps'
plot 'flint-factor/flint-timings.dat' using 1:7 with lines t 'Flint', \
     'pari-factor/pari-timings.dat' using 1:7 with lines t 'Pari', \
     'squfof-parigp/squfof-timings.dat' using 1:7 with lines t 'SQUFOF', \
     'spar-factor/spar-vanilla-timings.dat' using 1:7 with lines t 'Vanilla SPAR', \
     'sspar-factor/sspar-optimized-timings.dat' using 1:7 with lines t 'SuperSPAR'

# median timings
set xrange [*:*]
set output 'factor-median.eps'
plot 'flint-factor/flint-timings.dat' using 1:4 with lines t 'Flint', \
     'pari-factor/pari-timings.dat' using 1:4 with lines t 'Pari', \
     'squfof-parigp/squfof-timings.dat' using 1:4 with lines t 'SQUFOF', \
     'spar-factor/spar-vanilla-timings.dat' using 1:4 with lines t 'Vanilla SPAR', \
     'sspar-factor/sspar-optimized-timings.dat' using 1:4 with lines t 'SuperSPAR'

set xrange [44:58]
set output 'factor-median-zoom-left.eps'
plot 'flint-factor/flint-timings.dat' using 1:4 with lines t 'Flint', \
     'pari-factor/pari-timings.dat' using 1:4 with lines t 'Pari', \
     'squfof-parigp/squfof-timings.dat' using 1:4 with lines t 'SQUFOF', \
     'spar-factor/spar-vanilla-timings.dat' using 1:4 with lines t 'Vanilla SPAR', \
     'sspar-factor/sspar-optimized-timings.dat' using 1:4 with lines t 'SuperSPAR'

set xrange [58:70]
set output 'factor-median-zoom-right.eps'
plot 'flint-factor/flint-timings.dat' using 1:4 with lines t 'Flint', \
     'pari-factor/pari-timings.dat' using 1:4 with lines t 'Pari', \
     'squfof-parigp/squfof-timings.dat' using 1:4 with lines t 'SQUFOF', \
     'spar-factor/spar-vanilla-timings.dat' using 1:4 with lines t 'Vanilla SPAR', \
     'sspar-factor/sspar-optimized-timings.dat' using 1:4 with lines t 'SuperSPAR'