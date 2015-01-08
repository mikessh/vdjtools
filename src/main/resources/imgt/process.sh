IMGTPARSER="java -cp ../../../../target/vdjtools-1.0-SNAPSHOT.jar com.antigenomics.vdjtools.imgt.ImgtParser"
$IMGTPARSER -b F+ORF+in-frame_P_gaps_05-06-14.imgt segments
$IMGTPARSER -nb F+ORF+in-frame_P_gaps_05-06-14.imgt segments.all
$IMGTPARSER -nmb F+ORF+in-frame_P_gaps_05-06-14.imgt segments.all.minor