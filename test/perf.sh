PLACE="$HOME/sparse_large/matrices/"
rm -rf a stat

MATRICES="webbase-1M-symmetric.mtx dblp-2010.mtx cnr-2000-symmetric.mtx hollywood-2009.mtx"

for i in $MATRICES
do
    echo working on pagerank.jl $PLACE$i
    julia pagerank.jl $PLACE$i >& a
    echo  !!!! pagerank.jl $PLACE$i >> stat
    grep "seconds$\|perf$" a >> stat
done

MATRICES=" \
    audikw_1.mtx G3_circuit.mtx parabolic_fem.mtx thermal2.mtx tmt_sym_amd.mtx \
    inline_1.mtx thermomech_dM.mtx \
"

for i in $MATRICES
do
    echo working on pcg.jl $PLACE$i
    julia pcg.jl $PLACE$i >& a 
    echo  !!!! pcg.jl $PLACE$i >> stat 
    grep "seconds$\|perf$" a >> stat
done
