import os

def get_features():
    phenotype_file = config.get("phenotype_file")
    if not os.path.exists(phenotype_file):
        raise FileNotFoundError(f"Phenotype file {phenotype_file} not found.")
    with open(phenotype_file) as f:
        lines = f.readlines()[1:]  # Skip header
        feature_dict = {}
        for line in lines:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                print(f"Skipping malformed line: {line.strip()}")
                continue
            chr_ = parts[0].replace("chr", "")
            start = max(0, int(parts[1]) - 1000000) # This is the start of the cis window - not the start of the feature
            end = int(parts[2]) + 1000000 # This is the end of the cis window - not the end of the feature
            feature = parts[3].strip()
            feature_dict[feature] = {"chr": chr_, "start": start, "end": end}
    return feature_dict

feature_dict = get_features()
feature_list = list(feature_dict.keys())  # Extracting feature names (keys)

print("Features to process:", list(feature_dict.keys()))

rule all:
    input:
        expand("results/{feature}_clump.clumped", feature=feature_dict.keys())

rule extract_snps:
    input:
        bim_file=config["bim_file"]
    output:
        "results/{feature}_cis_region.snps"
    params:
        chr=lambda wildcards: feature_dict.get(wildcards.feature, {}).get("chr", "NA"),
        start=lambda wildcards: feature_dict.get(wildcards.feature, {}).get("start", 0),
        end=lambda wildcards: feature_dict.get(wildcards.feature, {}).get("end", 0)
    shell:
        """
        awk -v chr="{params.chr}" -v start={params.start} -v end={params.end} '
            $1 == chr && $4 >= start && $4 <= end {{ print $2 }}
        ' {input.bim_file} > {output}
        """

rule extract_sumstat:
    input:
        sumstat=config["sumstat_file"],  # Tabixed sumstat file
    output:
        "results/{feature}_sumstats.txt"
    params:
        chr=lambda wildcards: "chr" + str(feature_dict.get(wildcards.feature, {}).get("chr", "NA")).lstrip("chr"),
        start=lambda wildcards: feature_dict.get(wildcards.feature, {}).get("start", 0),
        end=lambda wildcards: feature_dict.get(wildcards.feature, {}).get("end", 0),
        columns=",".join(str(c) for c in config["columns_to_extract"]),
        col_names="\\t".join(config["column_names"])
    shell:
        """
        echo "Extracting sumstats for feature: {wildcards.feature}";
        mkdir -p $(dirname {output});

        ml tabix;

        # Extract specified columns and filter by feature
        tabix {input.sumstat} {params.chr}:{params.start}-{params.end} | \
        awk -F'\\t' -v feat='{wildcards.feature}' -v cols='{params.columns}' '
            BEGIN {{
                split(cols, col_nums, ",");
                print "{params.col_names}"
            }}
            $1 == feat {{
		for (i = 1; i <= length(col_nums); i++) {{
                    printf "%s%s", $(col_nums[i]+0), (i == length(col_nums) ? "\\n" : "\\t")
                }}
            }}
        ' > {output}

        # If the output is empty, add header only
        if [ $(wc -l < {output}) -le 1 ]; then
            echo -e "{params.col_names}" > {output}
        fi
        """

rule run_plink_clumping:
    input:
        sumstat="results/{feature}_sumstats.txt",
        region="results/{feature}_cis_region.snps"
    output:
        "results/{feature}_clump.clumped"
    params:
        plink_file=config["plink_file"],
        clump_p1=config["clump_p1"],
        clump_p2=config["clump_p2"],
        clump_r2=config["clump_r2"],
        clump_field=lambda wildcards: config["column_names"][3]  # Get the fourth column name
    shell:
        """
        ml plink/1.90
        plink --bfile {params.plink_file} \
              --clump {input.sumstat} \
              --maf 0.01 \
              --clump-snp-field SNP \
              --clump-field {params.clump_field} \
              --clump-p1 {params.clump_p1} \
              --clump-p2 {params.clump_p2} \
              --clump-r2 {params.clump_r2} \
              --extract {input.region} \
              --out results/{wildcards.feature}_clump || true
        
        # Check if PLINK output exists, if not create an empty file
        if [ ! -f "results/{wildcards.feature}_clump.clumped" ]; then
            echo -e "CHR\tSNP\tBP\tP" > results/{wildcards.feature}_clump.clumped
            echo "Warning: No significant --clump results. Empty file created."
        fi
        """

rule summarize_results:
    input:
        sumstats=expand("results/{feature}_sumstats.txt", feature=feature_list),
        clumped=expand("results/{feature}_clump.clumped", feature=feature_list),
        phenotype_file=config["phenotype_file"]  # Assuming feature info comes from here
    output:
        summary="results/snp_count_per_feature_clump_summary.tsv"
    shell:
        """
        echo -e "CHR\\tSTART\\tEND\\tFEATURE\\tCIS_REGION_SNP_COUNT\\tSNP_COUNT_CLUMPED" > {output.summary}

        while IFS=$'\t' read -r chr start end feature; do
            sumstat_file="results/${{feature}}_sumstats.txt"
            clumped_file="results/${{feature}}_clump.clumped"

            # Count SNPs in the sumstat file (excluding the header)
            if [ -f "$sumstat_file" ]; then
                snp_count=$(awk 'NR > 1 { count++ } END {{ print count ? count : 0 }}' "$sumstat_file")
            else
                snp_count=0
            fi

            # Count SNPs in the clumped file (excluding the header)
            if [ -f "$clumped_file" ]; then
                clumped_snp_count=$(awk 'NR > 1 { count++ } END {{ print count ? count : 0 }}' "$clumped_file")
            else
                clumped_snp_count=0
            fi

            # Append the results to the summary file
            echo -e "$chr\\t$start\\t$end\\t$feature\\t$snp_count\\t$clumped_snp_count" >> {output.summary}

        done < {input.phenotype_file}
        """

