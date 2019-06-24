set term latex
set output "q.tex"
factor=0.7
set size factor*5/5.,factor*4/3.
set format xy "$%g$"
set title  'Function $\mathsf{Q}(X)$'
set xlabel '$X$'
set ylabel '$\mathsf{Q}(X)$'
set nokey
set ytics 0,0.1,0.35
#set grid
plot [0:3] [0:0.4] x/(1+x+x**2)
