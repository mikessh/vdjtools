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
$VDJTOOLS PoolSamples -m metadata.small.txt out/13

# annotation
$VDJTOOLS Annotate -m metadata.txt out/annot/

# check all output files are generated

cd out/
ls -lh

flist=(
        '0.basicstats.txt'
        '1.spectratype.insert.wt.txt'
        '1.spectratype.ndn.wt.txt'
        '1.spectratype.nt.wt.txt'
        '2.segments.wt.J.pdf'
        '2.segments.wt.J.txt'
        '2.segments.wt.V.pdf'
        '2.segments.wt.V.txt'
        '3.fancyspectra.pdf'
        '3.fancyspectra.txt'
        '4.spectraV.wt.pdf'
        '4.spectraV.wt.txt'
        '5.fancyvj.wt.pdf'
        '5.fancyvj.wt.txt'
        '6.qstat.pdf'
        '6.qstat.txt'
        '7.diversity.strict.exact.txt'
        '7.diversity.strict.resampled.txt'
        '8.rarefaction.strict.pdf'
        '8.rarefaction.strict.txt'
        '9.paired.strict.summary.txt'
        '9.paired.strict.table.collapsed.pdf'
        '9.paired.strict.table.collapsed.txt'
        '9.paired.strict.table.txt'
        '9.strict.paired.scatter.pdf'
        '10.age.hc.aa.F.pdf'
        '10.age.mds.aa.F.pdf'
        '10.age.mds.aa.F.txt'
        '10.intersect.batch.aa.txt'
        '12.join.aa.summary.txt'
        '12.join.aa.table.txt'
        '12.join.aa.venn.pdf'
        '13.pool.aa.table.txt'
        'dec/metadata.txt'
        'ds/metadata.txt'
        'nf/metadata.txt'
        'annot/metadata.txt'
    )

for f in "${flist[@]}"
do
    if [[ ! -s $f ]]
        then exit 1
    fi
done