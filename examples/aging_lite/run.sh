VDJTOOLS="java -Xmx6G -jar ../../vdjtools-1.0-SNAPSHOT.jar"
PARAMS="-S mitcr -m samples/metadata.txt"

# Basic
$VDJTOOLS CalcBasicStats $PARAMS ./out/0
$VDJTOOLS CalcSpectratype $PARAMS ./out/1
$VDJTOOLS CalcSegmentUsage $PARAMS -p -f age -n ./out/2
$VDJTOOLS PlotFancySpectratype -S mitcr ./samples/A4-i125.txt.gz ./out/3
$VDJTOOLS PlotSpectratypeV -S mitcr ./samples/A4-i125.txt.gz ./out/4
$VDJTOOLS PlotFancyVJUsage -S mitcr ./samples/A4-i125.txt.gz ./out/5

# Diversity
$VDJTOOLS CalcDiversityStats $PARAMS ./out/6
$VDJTOOLS RarefactionPlot $PARAMS -f age -n -l sample.id ./out/7

# Intersect
$VDJTOOLS IntersectPair -S mitcr -c 100 -p ./samples/A4-i189.txt.gz ./samples/A4-i190.txt.gz ./out/8
$VDJTOOLS BatchIntersectPair $PARAMS ./out/9
$VDJTOOLS BatchIntersectPairPlot -f age -n -l sample.id ./out/9

# Annotation
$VDJTOOLS ScanDatabase $PARAMS -f --filter "__origin__=~/EBV/" ./out/10