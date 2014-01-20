set multiplot layout 2,1
set style line 1 lt 1 lc 1 lw 1 pt 6
set style line 2 lt 1 lc 3 lw 1 pt 6
set style line 3 lt 1 lc 5 lw 1 pt 6

set xlabel 'Simulation step'
set ylabel 'VdW energy [kcal/mol]'
set log x
set xrange [1:3000]
plot 'qsurr_benchmark.en' using 1:2                                                     notitle  with lines ls 1,            'qsurr_benchmark.en' using 1:3                                     title 'Benchmark bounds' with lines ls 1,            'this_qsurr.en' using 1:2                                          title 'The present build' with lines ls 2

set ylabel 'El energy [kcal/mol]'
plot 'qsurr_benchmark.en' using 1:4                                                          notitle  with lines ls 1,       'qsurr_benchmark.en' using 1:5                                          title 'Benchmark bounds' with lines ls 1,       'this_qsurr.en' using 1:3                                        title 'The present build' with lines ls 2     
