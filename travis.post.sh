VDJTOOLS="java -Xmx4G -jar ../target/vdjtools-*.jar"
cd aging_lite/

# basic analysis
$VDJTOOLS CalcBasicStats -m metadata.txt out/0
$VDJTOOLS CalcSpectratype -m metadata.txt out/1
$VDJTOOLS CalcSegmentUsage -m metadata.txt -p -f age -n out/2
$VDJTOOLS PlotFancySpectratype A4-i125.txt.gz out/3
$VDJTOOLS PlotSpectratypeV A4-i125.txt.gz out/4
$VDJTOOLS PlotFancyVJUsage A4-i125.txt.gz out/5

# diversity estimates
$VDJTOOLS PlotQuantileStats A4-i125.txt.gz out/6
$VDJTOOLS CalcDiversityStats -m metadata.txt out/7
$VDJTOOLS RarefactionPlot -m metadata.txt -f age -n -l sample.id out/8

# sample overlap
$VDJTOOLS OverlapPair -p A4-i189.txt.gz A4-i190.txt.gz out/9
$VDJTOOLS CalcPairwiseDistances -m metadata.small.txt out/10
$VDJTOOLS ClusterSamples -p -f age -n -l sample.id out/10 out/10.age

# sample operations and filtering
$VDJTOOLS Decontaminate -m metadata.txt -c out/dec/
$VDJTOOLS Downsample -m metadata.txt -c -x 10000 out/ds/
$VDJTOOLS FilterNonFunctional -m metadata.txt -c out/nf/
$VDJTOOLS JoinSamples -p -m metadata.small.txt out/12
$VDJTOOLS PoolSamples -w -m metadata.small.txt out/13

cd out/
ls -lh