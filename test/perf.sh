PLACE="~/sparse_large/matrices/"

MATRICES=" \
    audikw_1.mtx          G3_circuit.mtx      parabolic_fem.mtx  thermal2.mtx       tmt_sym_amd.mtx \
    delaunay_n24.mtx      hugetric-00020.mtx  rajat31.mtx        thermomech_dK.mtx  webbase-1M.mtx \
    dielFilterV3real.mtx  inline_1.mtx        road_usa.mtx       thermomech_dM.mtx
"

rm -rf a stat
for i in $MATRICES
do
    echo working on pcg.jl $PLACE$i
    julia pcg.jl $PLACE$i >& a 
    grep "seconds$\|original\|accelerated" a >> stat
done

for i in $MATRICES
do
    echo working on pagerank.jl $PLACE$i
	julia pagerank.jl $PLACE$i >& a
    grep "seconds$\|original\|accelerated" a >> stat
done

