import os
import subprocess
import logging
import argparse
import pysam
import csv
import glob
import shutil

def extract_junctions(bam_list_path, junction_output_dir, junction_file_list):
    """
    Reads a TSV with:
      col1 = BAM path (required)
      col2 = group (optional but recommended for DS)
    Writes:
      - <output_dir>/junction_files/*.junc
      - <output_dir>/junction_files.txt   (newline-separated .junc paths)
      - <output_dir>/samples.tsv          (sample or sample\tgroup)
    Returns:
      (junction_list_file, samples_tsv)
    """
    os.makedirs(junction_output_dir, exist_ok=True)

    out_dir = os.path.dirname(junction_file_list)
    samples_tsv = os.path.join(out_dir, "samples.tsv")
    group_file_tsv = os.path.join(out_dir, "group_file.tsv")

    with open(bam_list_path, "r") as f_in, open(junction_file_list, "w") as f_junc_list, open(samples_tsv, "w") as f_samples:

        reader = csv.reader(f_in, delimiter="\t")

        for row in reader:
            if not row or not row[0].strip() or row[0].lstrip().startswith("#"):
                continue

            bam_file = row[0].strip()
            sample_id = os.path.basename(bam_file).replace(".bam", "")
            group = row[1].strip() if len(row) > 1 and row[1].strip() else None

            # Write samples.tsv (1 or 2 columns)
            if group:
                f_samples.write(f"{sample_id}\t{group}\n")
            else:
                f_samples.write(f"{sample_id}\n")

            # # Always write a 2-col group file (needed by leafcutter_ds.R)
            # f_group.write(f"{sample_id}\t{group if group else 'UNKNOWN'}\n")

            junc_file = os.path.join(junction_output_dir, f"{sample_id}.junc")

            if os.path.exists(junc_file):
                logging.info(f"Junction file exists; skipping: {junc_file}")
                f_junc_list.write(junc_file + "\n")
                continue

            try:
                # creates .bai alongside the BAM; will fail if BAM location is read-only
                pysam.index(bam_file)
            except Exception as e:
                logging.error(f"Failed to index BAM file {bam_file}: {e}")
                continue

            logging.info(f"Extracting junctions from {bam_file} -> {junc_file}")
            cmd = [
                "regtools", "junctions", "extract",
                "-a", "8",
                "-m", "50",
                "-M", "500000",
                "-s", "FR",
                bam_file,
                "-o", junc_file
            ]
            subprocess.run(cmd, check=True)
            f_junc_list.write(junc_file + "\n")

    logging.info(f"Wrote: {junction_file_list}, {samples_tsv}")
    return junction_file_list, samples_tsv


def cluster_introns(junctions_file_list, out_dir, cell_type,
                    script="leafcutter/clustering/leafcutter_cluster_regtools.py",
                    py="python"):
    os.makedirs(out_dir, exist_ok=True)

    # Leafcutter writes to CWD using prefix from -o
    cmd = [
        py, script,
        "-j", junctions_file_list,
        "-o", cell_type,
        "-m", "50",
        "-l", "500000"
    ]
    subprocess.run(cmd, check=True, cwd=out_dir)

    # Tidy per-sample cluster files
    per_sample_clusters_dir = os.path.join(out_dir, "per_sample_clusters")
    os.makedirs(per_sample_clusters_dir, exist_ok=True)
    for file in glob.glob(os.path.join(out_dir, f"*.junc.{cell_type}.sorted.gz")):
        shutil.move(file, os.path.join(per_sample_clusters_dir, os.path.basename(file)))

    logging.info(f"Intron clustering finished. Prefix: {cell_type}")


def differential_quantifications(counts_file, group_file, out_dir, exon_file,
                                 min_samples_per_intron, cell_type,
                                 rscript="Rscript",
                                 script_path="leafcutter/scripts/leafcutter_ds.R"):
    os.makedirs(out_dir, exist_ok=True)
    output_prefix = os.path.join(out_dir, cell_type)

    cmd = [
        rscript, script_path,
        counts_file,
        group_file,
        "-o", output_prefix,
        "-i", str(min_samples_per_intron),
        "-e", exon_file
    ]
    subprocess.run(cmd, check=True, cwd=out_dir)
    logging.info(f"Differential splicing finished. Results at prefix: {output_prefix}")


def main():
    parser = argparse.ArgumentParser(description="LeafCutter pipeline")
    parser.add_argument('--bam_file_list', type=str, required=True, help="TSV: BAM path [TAB group]")
    parser.add_argument('--output_dir', type=str, required=True, help="Output directory")
    parser.add_argument('--cell_type', type=str, required=True, help="Cell type label")
    parser.add_argument('--exon_file', type=str, required=True, help="Exon file (BED/GTF-derived) for DS")
    args = parser.parse_args()

    # ---- logging FIRST  ----
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)
    log_dir  = os.path.join(output_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, f"{args.cell_type}.log")
    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s: %(message)s',
        handlers=[logging.FileHandler(log_file, mode='w'), logging.StreamHandler()],
        force=True,
    )
    logging.info(f"Logging to {log_file}")

    # ---- paths ----
    junction_output_dir = os.path.join(output_dir, "junction_files")
    junction_list_file  = os.path.join(output_dir, "junction_files.txt")

    # ---- Step 1: Extract junctions + TSVs ----
    # Make sure extract_junctions returns (jlist, samples_tsv, group_file_tsv)
    jlist, samples_tsv = extract_junctions(
        args.bam_file_list, junction_output_dir, junction_list_file
    )

    # ---- Step 2: Clustering ----
    cluster_introns(jlist, output_dir, args.cell_type)

    # Leafcutter writes <prefix>_perind_numers.counts.gz into out_dir
    counts_file = os.path.join(output_dir, f"{args.cell_type}_perind_numers.counts.gz")
    if not os.path.isfile(counts_file):
        logging.error(f"Counts file not found: {counts_file}")
        raise FileNotFoundError(counts_file)

    # ---- Step 3: Differential splicing ----
    # Use the 2-col group file, NOT samples.tsv
    differential_quantifications(
        counts_file, samples_tsv, output_dir, args.exon_file, 3, args.cell_type
    )

if __name__ == "__main__":
    main()

# if __name__ == "__main__":
#     logging.basicConfig(
#         format='[%(asctime)s] %(levelname)s: %(message)s',
#         level=logging.INFO
#     )
#     main()
