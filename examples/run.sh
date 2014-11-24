VDJTOOLS="java -Xmx6G -jar ../target/vdjtools-1.0-SNAPSHOT.jar"
PARAMS="-S mitcr -m samples/metadata.txt"

# Basic
#$VDJTOOLS CalcBasicStats $PARAMS ./out/0
#$VDJTOOLS CalcSpectratype $PARAMS ./out/1
#$VDJTOOLS PlotFancySpectratype -S mitcr ./samples/ds10k.A4-i127.txt.gz ./out/2
#$VDJTOOLS PlotSpectratypeV -S mitcr ./samples/ds10k.A4-i127.txt.gz ./out/3
#$VDJTOOLS CalcSegmentUsage $PARAMS -p -f age -n ./out/4
#$VDJTOOLS PlotFancyVJUsage -S mitcr ./samples/ds10k.A4-i127.txt.gz ./out/5

# Sample size and diversity
#$VDJTOOLS CalcDiversityStats -S mitcr -m metadata_small.txt -n 5000 ./out/6
#$VDJTOOLS RarefactionPlot -S mitcr -m metadata_small.txt -r 3 -s 10 -l sex -f age -n ./out/6
#$VDJTOOLS CalcDiversityStats -S mitcr -m metadata.txt -r -r-max 10000 -r-step 500 -n 5000 -p ./out/6
#$VDJTOOLS DownSample -S mitcr -m ds/metadata.txt -n 100 ./ds/2/2
#$VDJTOOLS BuildFrequencyTable -p -S mitcr ./samples/ds10k.A5-S1.txt.gz ./out/8

# Intersect
#$VDJTOOLS IntersectPair -S mitcr -c 100 -p ./sample/ds10k.A6-I208ob.txt.gz ./sample/ds10k.A5-S6.txt.gz ./out/9
#$VDJTOOLS BatchIntersectPair $PARAMS ./out/10
#$VDJTOOLS BatchIntersectPairPlot -m vJSD -f sex ./out/10.batch_intersect_aa.txt ./out/10
#$VDJTOOLS IntersectSequential -S simple -m metadata_tc.txt -x 0 -c 100 -p ./out/11
#$VDJTOOLS PoolSamples $PARAMS -w ./out/12


# Annot
$VDJTOOLS ScanDatabase $PARAMS -f --filter "__pathogen__=~/EBV/" ./out/13