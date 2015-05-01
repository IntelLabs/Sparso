# Note: the directory in PLACE must be an absolute path. Otherwise, segmentation fault
PLACE="/home/jpark103/matrices/"
rm -rf a stat

# this is the setting for endevor hsw arch. change OMP_NUM_THREADS to the number of the cores of your own machine
export OMP_NUM_THREADS=14
export KMP_AFFINITY=granularity=fine,compact,1
 
 
#MATRICES="webbase-1M-symmetric.mtx dblp-2010.mtx cnr-2000-symmetric.mtx hollywood-2009.mtx"
#for i in $MATRICES
 
while read i; 
do
    echo working on pagerank.jl $PLACE$i.mtx
    julia pagerank.jl $PLACE$i.mtx
    echo
done < pagerank.lst
  
while read i;
do
    echo working on pcg.jl $PLACE$i.mtx
    julia pcg.jl $PLACE$i.mtx
    echo
done < pcg.lst
