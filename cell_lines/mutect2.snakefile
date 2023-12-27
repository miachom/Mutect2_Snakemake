configfile: "config/samples.yaml"
configfile: "config/config.yaml"

rule all:
        input:
                expand("results/LearnReadOrientationModel/{tumor}/read_orientation_model.tar.gz",tumor=config["normals"]),
                expand("results/GetPileupSummaries/{tumor}/pileup_summaries.table",tumor=config["normals"]),
                expand("results/CalculateContamination/{tumor}/{tumor}_contamination.table",tumor=config["normals"]),
                expand("results/CalculateContamination/{tumor}/{tumor}.segments.table",tumor=config["normals"]),
                expand("results/FilterMutectCalls/{tumor}/filtered_all.vcf.gz",tumor=config["normals"]),
                expand("results/FilterMutectCalls/{tumor}/filtering_stats.tsv",tumor=config["normals"])

rule mutect2:
        input:
                tumor_filepath = lambda wildcards: config["base_file_name"][wildcards.tumor]
        output:
                vcf = "results/mutect2/{tumor}/{tumor}_unfiltered.vcf.gz",
                tbi = "results/mutect2/{tumor}/{tumor}_unfiltered.vcf.gz.tbi",
                tar = "results/mutect2/{tumor}/{tumor}_unfiltered.f1r2.tar.gz",
                stats = "results/mutect2/{tumor}/{tumor}_unfiltered.vcf.gz.stats",
                bam = "results/mutect2/{tumor}/{tumor}_unfiltered.bamout.bam",
                bai = "results/mutect2/{tumor}/{tumor}_unfiltered.bamout.bai"
        params:
                reference_genome = config["reference_genome"],
                germline_resource = config["germline_resource"],
                gatk = config["gatk_path"],
                panel_of_normals = config["panel_of_normals"],
                normals = lambda wildcards: config["normals"][wildcards.tumor],
                intervals = config["interval_list"]
        log:
                "logs/mutect2/{tumor}_mutect2.txt"
        shell:
                "({params.gatk} Mutect2 \
                -reference {params.reference_genome} \
                -input {input.tumor_filepath} \
                -tumor {params.normals} \
                -intervals {params.intervals} \
                --interval-padding 100 \
                --germline-resource {params.germline_resource} \
                --genotype-germline-sites true \
                --genotype-pon-sites true \
                --f1r2-tar-gz {output.tar} \
                --panel-of-normals {params.panel_of_normals} \
                --read-filter PassesVendorQualityCheckReadFilter \
                --read-filter HasReadGroupReadFilter \
                --read-filter NotDuplicateReadFilter \
                --read-filter MappingQualityAvailableReadFilter \
                --read-filter MappingQualityReadFilter \
                --minimum-mapping-quality 30 \
                --read-filter OverclippedReadFilter \
                --filter-too-short 25 \
                --read-filter GoodCigarReadFilter \
                --read-filter AmbiguousBaseReadFilter \
                --native-pair-hmm-threads 2 \
                --seconds-between-progress-updates 100 \
                --downsampling-stride 20 \
                --max-reads-per-alignment-start 6 \
                --max-suspicious-reads-per-alignment-start 6 \
                --bam-output {output.bam} \
                -output {output.vcf}) 2> {log}"

rule LearnReadOrientationModel:
        input:
                "results/mutect2/{tumor}/{tumor}_unfiltered.f1r2.tar.gz"
        output:
                "results/LearnReadOrientationModel/{tumor}/read_orientation_model.tar.gz"
        params:
                gatk = config["gatk_path"]
        log:
                "logs/LearnReadOrientationModel/{tumor}/LearnReadOrientationModel.txt"
        shell:
                "({params.gatk} LearnReadOrientationModel \
                -I {input} \
                -O {output}) 2> {log}"

rule GetPileupSummaries:
        input:
                filepath = lambda wildcards: config["base_file_name"][wildcards.tumor]
                output:
                "results/GetPileupSummaries/{tumor}/pileup_summaries.table"
        params:
                reference_genome = config["reference_genome"],
                gatk = config["gatk_path"],
                variants_for_contamination = config["variants_for_contamination"],
                intervals = config["interval_list"]
        log:
                "logs/GetPileupSummaries/{tumor}/GetPileupSummaries.txt"
        shell:
                """
                ({params.gatk} GetPileupSummaries \
                -R {params.reference_genome} \
                -I {input.filepath} \
                -V {params.variants_for_contamination} \
                -L {params.intervals} \
                -O {output}) 2> {log}
                """

rule CalculateContamination:
input:
                tumor_pileup="results/GetPileupSummaries/{tumor}/pileup_summaries.table"
        output:
                contamination_table="results/CalculateContamination/{tumor}/{tumor}_contamination.table",
                tumor_segmentation="results/CalculateContamination/{tumor}/{tumor}.segments.table"
        params:
                gatk = config["gatk_path"]
        log:
                "logs/CalculateContamination/{tumor}/contamination.log"
        shell:
                "({params.gatk} CalculateContamination \
                -I {input.tumor_pileup} \
                --tumor-segmentation {output.tumor_segmentation} \
                -O {output.contamination_table}) 2> {log}"

rule FilterMutectCalls:
        input:
                unfiltered_vcf = "results/mutect2/{tumor}/{tumor}_unfiltered.vcf.gz",
                vcf_index = "results/mutect2/{tumor}/{tumor}_unfiltered.vcf.gz.tbi",
                segments_table = "results/CalculateContamination/{tumor}/{tumor}.segments.table",
                contamination_table = "results/CalculateContamination/{tumor}/{tumor}_contamination.table",
                read_orientation_model = "results/LearnReadOrientationModel/{tumor}/read_orientation_model.tar.gz",
                mutect_stats = "results/mutect2/{tumor}/{tumor}_unfiltered.vcf.gz.stats"
        output:
                filtered_vcf = "results/FilterMutectCalls/{tumor}/filtered_all.vcf.gz",
                filtering_stats = "results/FilterMutectCalls/{tumor}/filtering_stats.tsv"
        params:
                gatk = config["gatk_path"],
                reference_genome = config["reference_genome"]
        log:
                "logs/FilterMutectCalls/{tumor}/FilterMutectCalls"
        shell:
                "({params.gatk} FilterMutectCalls \
                -R {params.reference_genome} \
                -V {input.unfiltered_vcf} \
                --tumor-segmentation {input.segments_table} \
                --contamination-table {input.contamination_table} \
                --ob-priors {input.read_orientation_model} \
                --stats {input.mutect_stats} \
                --filtering-stats {output.filtering_stats} \
                -O {output.filtered_vcf}) 2> {log}"
