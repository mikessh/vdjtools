VDJTOOLS="java -Xmx6G -jar ../../target/vdjtools-1.0-SNAPSHOT.jar"
PARAMS="-S simple -m samples/metadata.txt"

# Basic
$VDJTOOLS CalcBasicStats $PARAMS ./out/0
$VDJTOOLS CalcSpectratype $PARAMS ./out/1
$VDJTOOLS CalcSegmentUsage $PARAMS -p -f "Time post HSCT, months" -n ./out/2

# Diversity
$VDJTOOLS CalcDiversityStats $PARAMS ./out/3
$VDJTOOLS RarefactionPlot $PARAMS -f "Time post HSCT, months" -n -l sample.id ./out/4

# Intersect
$VDJTOOLS IntersectPair -S simple -p ./samples/minus48months.txt.gz ./samples/4months.txt.gz ./out/5
$VDJTOOLS IntersectSequential $PARAMS -f "Time post HSCT, months" -x 0 -p ./out/6

# Annotation
$VDJTOOLS ScanDatabase $PARAMS -f --filter "__origin__.contains('CMV')||__origin__.contains('EBV')" ./out/7

# Manipulation
$VDJTOOLS ApplySampleAsFilter $PARAMS -c ./samples/minus48months.txt.gz out/asaf/