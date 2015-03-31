PLACE="$HOME/sparse_large/matrices/"
#PLACE="../../../sparse_large/matrices/"

while read p; do
    echo CSR_test $p ....
    ./CSR_test $PLACE/$p.mtx 
done < ../../test/pagerank.lst

while read p; do
    echo CSR_test $p ....
    ./CSR_test $PLACE/$p.mtx 
done < ../../test/pcg.lst
