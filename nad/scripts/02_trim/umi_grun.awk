# NAD branch step 2, second half: runs on the FASTQ that cutadapt writes after
# removing the 3' adapter, and does the three things cutadapt cannot do in one
# pass.
#
# Read structure (GSE146368 NAD captureSeq, checked on SRR11236402/444):
#   [G run from the TdT C-tail] [RNA insert] [C] [NNNNNN UMI] [AGATCGGAAGAGC...]
#
# cutadapt is run with --rename '{id} ad={adapter_name}', so the header tells
# whether the adapter was found. For each read this script:
#   1. if the adapter was found, removes the last 7 bases (C + 6-nt UMI) and
#      keeps the UMI; otherwise the UMI is unknown and recorded as NNNNNN;
#   2. strips the leading G run and records its length;
#   3. drops reads whose remaining insert is shorter than MINLEN.
# The new read name is  <id>_g<Glen>_<UMI>  with no spaces, because STAR keeps
# only the first word of a name and umi_tools takes the UMI after the last '_'.
#
# Per-read G-run lengths are tallied and written to GRUN_TSV as
#   adapter  g_run  reads
# so the leading-G distribution can be inspected per library.
#
# Variables: MINLEN (default 18), GRUN_TSV (required), STATS_TSV (required).

BEGIN {
    if (MINLEN == "") MINLEN = 18
    if (GRUN_TSV == "" || STATS_TSV == "") { print "GRUN_TSV and STATS_TSV must be set" > "/dev/stderr"; exit 2 }
}

NR % 4 == 1 { hdr = $0; next }
NR % 4 == 2 { seq = $0; next }
NR % 4 == 3 { next }
NR % 4 == 0 {
    qual = $0
    id = substr(hdr, 2); sub(/[ \t].*/, "", id)
    has_ad = (hdr ~ / ad=no_adapter/) ? 0 : 1
    n_in++

    umi = "NNNNNN"
    if (has_ad) {
        L = length(seq)
        if (L < 7) { n_short++; next }
        umi = substr(seq, L - 5, 6)
        seq = substr(seq, 1, L - 7); qual = substr(qual, 1, L - 7)
    }

    match(seq, /^G*/); g = RLENGTH
    seq = substr(seq, g + 1); qual = substr(qual, g + 1)
    grun[has_ad "\t" g]++

    if (length(seq) < MINLEN) { n_short++; next }
    n_out++; if (has_ad) n_out_ad++
    printf "@%s_g%d_%s\n%s\n+\n%s\n", id, g, umi, seq, qual
}

END {
    print "adapter\tg_run\treads" > GRUN_TSV
    for (k in grun) print k "\t" grun[k] > GRUN_TSV
    print "reads_in\treads_out\treads_out_with_umi\treads_too_short" > STATS_TSV
    printf "%d\t%d\t%d\t%d\n", n_in, n_out, n_out_ad, n_short > STATS_TSV
}
