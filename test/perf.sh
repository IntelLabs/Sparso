# Note: the directory in PLACE must be an absolute path. Otherwise, segmentation fault
PLACE="/home/jpark103/matrices/"
rm -rf a stat
 
#MATRICES="webbase-1M-symmetric.mtx dblp-2010.mtx cnr-2000-symmetric.mtx hollywood-2009.mtx"
#for i in $MATRICES
 
while read i; 
do
    echo working on pagerank.jl $PLACE$i.mtx
    julia pagerank.jl $PLACE$i.mtx >& a
    echo  !!!! pagerank.jl $PLACE$i.mtx >> stat
    grep "seconds$\|perf$" a >> stat
done < pagerank.lst
  
while read i;
do
    echo working on pcg.jl $PLACE$i.mtx
    julia pcg.jl $PLACE$i.mtx >& a
    echo  !!!! pcg.jl $PLACE$i.mtx >> stat
    grep "seconds$\|perf$" a >> stat
done < pcg.lst
