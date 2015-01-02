VDJTOOLS="java -Xmx6G -jar ../../target/vdjtools-1.0-SNAPSHOT.jar"
PARAMS="-S mitcr -m samples/metadata.txt"

# Basic
$VDJTOOLS CalcBasicStats $PARAMS ./out/0
$VDJTOOLS CalcSpectratype $PARAMS ./out/1
$VDJTOOLS CalcSegmentUsage $PARAMS -p -f age -n ./out/2
$VDJTOOLS PlotFancySpectratype -S mitcr ./samples/A4-i125.txt.gz ./out/3
$VDJTOOLS PlotSpectratypeV -S mitcr ./samples/A4-i125.txt.gz ./out/4
$VDJTOOLS PlotFancyVJUsage -S mitcr ./samples/A4-i125.txt.gz ./out/5

# Diversity
$VDJTOOLS PlotQuantileStats -S mitcr ./samples/A4-i125.txt.gz ./out/6
$VDJTOOLS CalcDiversityStats $PARAMS -x 5000 ./out/7
$VDJTOOLS RarefactionPlot $PARAMS -f age -n -l sample.id ./out/8

# Intersect
$VDJTOOLS IntersectPair -S mitcr -p ./samples/A4-i189.txt.gz ./samples/A4-i190.txt.gz ./out/9
$VDJTOOLS BatchIntersectPair $PARAMS ./out/10
$VDJTOOLS BatchIntersectPairPlot -f age -n -l sample.id ./out/10 ./out/10.age
$VDJTOOLS BatchIntersectPairPlot -m vJSD -f sex -l sample.id ./out/10 ./out/10.sex

# Annotation
$VDJTOOLS ScanDatabase $PARAMS -f --filter "__origin__=~/EBV/" ./out/11

# Manipulation
$VDJTOOLS Decontaminate $PARAMS -c out/dec/
$VDJTOOLS Downsample $PARAMS -c -x 5000 out/ds/
$VDJTOOLS FilterNonFunctional $PARAMS -c out/nf/