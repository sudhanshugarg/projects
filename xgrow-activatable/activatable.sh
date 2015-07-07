actualSize=32
array=`expr $actualSize \- 1`
maxSize=`expr $array \* $array`
kf=3500000
events=`expr $actualSize \* $actualSize \* 2 \* 60 \* 10000`
GmcMax=100
GseMax=1
Gse1=26.57
Gse3=31.41
Gse5=36.24
Gse7=41.07
Gse9=45.90

for mc in $(seq 1 $GmcMax)
do
    for se in 1 3 5 7 9
    do
        #tmc=`echo "$events*e($mc)/$kf" | bc -l`
        #echo "./xgrow tilesets/my.sierpinski.tiles Gmc=$mc Gse=$se datafile=otm.dat emax=$events fsmax=$maxSize size=$actualSize k=$kf -nw"
        var="Gse$se"
        echo "./xgrow tilesets/activatable.$se-1$se.sierpinski.tiles Gmc=$mc Gse=${!var} datafile=aKTAM-$se-1$se.dat emax=$events fsmax=$maxSize size=$actualSize k=$kf -nw"
        ./xgrow tilesets/activatable.$se-1$se.sierpinski.tiles Gmc=$mc Gse=${!var} datafile=aKTAM-$se-1$se.dat emax=$events fsmax=$maxSize size=$actualSize k=$kf -nw
        echo "./xgrow_orig tilesets/sierpinski.tiles Gmc=$mc Gse=${!var} datafile=kTAM-1$se.dat emax=$events fsmax=$maxSize size=$actualSize k=$kf -nw"
        ./xgrow_orig tilesets/sierpinski.tiles Gmc=$mc Gse=${!var} datafile=kTAM-1$se.dat emax=$events fsmax=$maxSize size=$actualSize k=$kf -nw
    done
done
