actualSize=32
array=`expr $actualSize \- 1`
maxSize=`expr $array \* $array`
kf=3500000
events=`expr $actualSize \* $actualSize \* 2 \* 60`
GmcMax=80
GseMax=1

for mc in $(seq 1 $GmcMax)
do
    for se in $(seq 1 $GseMax)
    do
        #tmc=`echo "$events*e($mc)/$kf" | bc -l`
        #echo "./xgrow tilesets/my.sierpinski.tiles Gmc=$mc Gse=$se datafile=otm.dat emax=$events fsmax=$maxSize size=$actualSize k=$kf -nw"
#        echo "./xgrow tilesets/my.sierpinski.tiles Gmc=$mc Gse=36.24 datafile=otm.dat emax=1000000 fsmax=$maxSize size=$actualSize k=$kf -nw"
#        ./xgrow tilesets/my.sierpinski.tiles Gmc=$mc Gse=36.24 datafile=otm.dat emax=1000000 fsmax=$maxSize size=$actualSize k=$kf -nw
        echo "./xgrow_orig tilesets/sierpinski.tiles Gmc=$mc Gse=36.24 datafile=sierpinski-compare.dat emax=1000000 fsmax=$maxSize size=$actualSize k=$kf -nw"
        ./xgrow_orig tilesets/sierpinski.tiles Gmc=$mc Gse=36.24 datafile=sierpinski-compare.dat emax=1000000 fsmax=$maxSize size=$actualSize k=$kf -nw
    done
done
