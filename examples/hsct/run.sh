VDJTOOLS="java -Xmx6G -jar ../../target/vdjtools-1.0-SNAPSHOT.jar"
PARAMS="-S simple -m samples/metadata.txt"

# Basic
$VDJTOOLS CalcBasicStats $PARAMS ./out/0
$VDJTOOLS CalcSpectratype $PARAMS ./out/1
$VDJTOOLS CalcSegmentUsage $PARAMS -p -f time -n ./out/2

# Diversity
$VDJTOOLS CalcDiversityStats $PARAMS ./out/3
$VDJTOOLS RarefactionPlot $PARAMS -f time -n -l sample.id ./out/4

# Intersect
$VDJTOOLS IntersectSequential $PARAMS -c 100 -p ./out/5

# Annotation
$VDJTOOLS ScanDatabase $PARAMS -f --filter "__origin__.contains('CMV')||__origin__.contains('EBV')" ./out/6