actualSize=32
array=`expr $actualSize \- 1`
maxSize=`expr $array \* $array`
kf=600000
events=`expr $actualSize \* $actualSize \* 2 \* 60`
GmcMax=30
GseMax=30

for mc in $(seq 1 $GmcMax)
do
    for se in $(seq 1 $GseMax)
    do
        #tmc=`echo "$events*e($mc)/$kf" | bc -l`
        echo "./xgrow tilesets/my.sierpinski.tiles Gmc=$mc Gse=$se datafile=otm.dat emax=$events fsmax=$maxSize size=$actualSize k=$kf -nw"
        ./xgrow tilesets/my.sierpinski.tiles Gmc=$mc Gse=$se datafile=otm.dat emax=$events fsmax=$maxSize size=$actualSize k=$kf -nw
    done
done
