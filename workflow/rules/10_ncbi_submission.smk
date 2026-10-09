# -----------------------------------------------------------------------------
# NCBI GenBank submission, for genomes with ncbi_submission set.
#
# funannotate2 annotate runs table2asn too, but with a placeholder submitter
# (a hard-coded .sbt), and without table2asn's genome mode, so its .gbk is
# not submittable. ncbi_submission takes annotate's .tbl and the unmasked
# assembly, puts them in the form NCBI wants (prepare_ncbi_submission.py),
# and runs NCBI's table2asn with:
#   -t        your submission template (ncbi.sbt_template)
#   -M n      genome mode: validation + discrepancy report with FATAL marks
#   -euk -T   eukaryote; look the organism up in NCBI Taxonomy (network),
#             without which table2asn applies bacterial checks
#   -j        organism, strain, nuclear genetic code, tech=wgs, BioProject,
#             BioSample; per-sequence chromosome/organelle via -src-file
#   -gaps-min/-l  runs of Ns become assembly_gap features
# check_ncbi_validation.py then fails the rule if the validator reports
# ERROR/REJECT messages or the discrepancy report has FATAL categories
# (unless ncbi.allow_errors). table2asn's own files stay in
# 08_ncbi/table2asn/out/ either way, for looking into what went wrong.
#
# Upload 08_ncbi/<genome>.sqn through the Genome Submission Portal
# (for "update": as an update of the existing genome).
# -----------------------------------------------------------------------------
rule get_table2asn:
    output:
        table2asn = protected(TABLE2ASN),
    params:
        url = config["ncbi"]["table2asn_url"],
        sha256 = config["ncbi"]["table2asn_sha256"],
        dir = os.path.dirname(TABLE2ASN),
    log:
        "logs/get_table2asn.log",
    resources:
        mem_mb = resources["get_table2asn"]["mem_mb"],
        runtime = resources["get_table2asn"]["runtime"],
    container:
        PYTHON_CONTAINER
    shell:
        r"""
        mkdir -p {params.dir}
        python workflow/scripts/download.py "{params.url}" {params.dir}/table2asn.gz \
            --sha256 {params.sha256} > {log} 2>&1
        gunzip -f {params.dir}/table2asn.gz
        chmod +x {output.table2asn}
        {output.table2asn} -version >> {log} 2>&1
        """


def _source_modifiers(genome):
    meta = config["genomes"][genome]
    mods = [f"[organism={meta['species']}]"]
    if meta.get("strain"):
        mods.append(f"[strain={meta['strain']}]")
    mods += ["[gcode=1]", "[tech=wgs]"]
    if meta.get("bioproject"):
        mods.append(f"[bioproject={meta['bioproject']}]")
    if meta.get("biosample"):
        mods.append(f"[biosample={meta['biosample']}]")
    return " ".join(mods)


rule ncbi_submission:
    input:
        table2asn = TABLE2ASN,
        # The unmasked assembly; annotate ran on the soft-masked copy, which
        # has the same sequences and coordinates.
        genome = rules.funannotate2_clean.output,
        tbl = rules.funannotate2_annotate.output.tbl,
        curate = rules.funannotate2_annotate.output.curate,
        seqinfo = rules.sequence_info.output,
    output:
        sqn = "results/{genome}/08_ncbi/{genome}.sqn",
        gbf = "results/{genome}/08_ncbi/{genome}.gbf",
        report = "results/{genome}/08_ncbi/{genome}.validation.txt",
    params:
        work = "results/{genome}/08_ncbi/table2asn",
        sbt = config["ncbi"]["sbt_template"],
        comment = config["ncbi"].get("structured_comment") or "",
        modifiers = lambda wc: _source_modifiers(wc.genome),
        mito_gcode = lambda wc: config["genomes"][wc.genome].get("mito_gcode", 4),
        gaps_min = config["ncbi"].get("gaps_min", 10),
        linkage = config["ncbi"].get("linkage_evidence", "paired-ends"),
        allow_errors = "--allow-errors" if config["ncbi"].get("allow_errors") else "",
        extra = config["ncbi"].get("extra", ""),
    log:
        "logs/{genome}/ncbi_submission.log",
    resources:
        mem_mb = resources["ncbi_submission"]["mem_mb"],
        runtime = resources["ncbi_submission"]["runtime"],
    container:
        PYTHON_CONTAINER
    shell:
        r"""
        if [ ! -s "{params.sbt}" ]; then
            echo "ERROR: no NCBI submission template at {params.sbt} (ncbi.sbt_template)." \
                 "Make one at https://submit.ncbi.nlm.nih.gov/genbank/template/submission/" | tee {log} >&2
            exit 1
        fi
        rm -rf {params.work}
        mkdir -p {params.work}

        python workflow/scripts/prepare_ncbi_submission.py \
            --fasta {input.genome} \
            --tbl {input.tbl} \
            --seqinfo {input.seqinfo} \
            --indir {params.work}/in \
            --name {wildcards.genome} \
            --src {params.work}/source.tsv \
            --mito-gcode {params.mito_gcode} \
            > {log} 2>&1

        src_flag=()
        [ -s {params.work}/source.tsv ] && src_flag=(-src-file {params.work}/source.tsv)
        comment_flag=()
        [ -n "{params.comment}" ] && comment_flag=(-w "{params.comment}")

        {input.table2asn} \
            -t {params.sbt} \
            -indir {params.work}/in \
            -outdir {params.work}/out \
            -M n -Z -euk -T -V b \
            -gaps-min {params.gaps_min} \
            -l {params.linkage} \
            -j "{params.modifiers}" \
            "${{src_flag[@]}}" \
            "${{comment_flag[@]}}" \
            {params.extra} \
            >> {log} 2>&1

        python workflow/scripts/check_ncbi_validation.py \
            --outdir {params.work}/out \
            --report {params.work}/validation.txt \
            --curate {input.curate} \
            {params.allow_errors} \
            >> {log} 2>&1

        cp {params.work}/out/{wildcards.genome}.sqn {output.sqn}
        cp {params.work}/out/{wildcards.genome}.gbf {output.gbf}
        cp {params.work}/validation.txt {output.report}
        """
