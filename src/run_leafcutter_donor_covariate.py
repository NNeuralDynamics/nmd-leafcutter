import os
import sys
import gzip
import subprocess
import logging
import argparse
import csv
import glob
import shutil

import pysam

# All default paths are relative to THIS script, so the pipeline works
# from any folder. Expected layout (either works):
#   nmd-leafcutter/run_leafcutter_donor_covariate.py + nmd-leafcutter/leafcutter/
#   nmd-leafcutter/src/run_leafcutter_donor_covariate.py + nmd-leafcutter/leafcutter/
HERE = os.path.dirname(os.path.abspath(__file__))


def default_leafcutter_dir():
    for cand in (os.path.join(HERE, "leafcutter"), os.path.join(HERE, "..", "leafcutter")):
        if os.path.isdir(cand):
            return os.path.normpath(cand)
    return os.path.join(HERE, "leafcutter")


def extract_junctions(bam_list_path, junction_output_dir, junction_file_list, strand="RF"):
    """
    Reads a headerless, tab-separated BAM list:
      col1 = BAM path (required; relative paths are read from the folder you launch in)
      col2 = treatment / group (required for DS)
      col3 = donor (optional, used as a covariate by leafcutter_ds.R)
    Writes:
      - <output_dir>/junction_files/*.junc
      - <output_dir>/junction_files.txt   (one .junc path per line)
      - <output_dir>/samples.tsv          (name, treatment, donor)
    Returns:
      (junction_list_file, samples_tsv)
    """
    os.makedirs(junction_output_dir, exist_ok=True)

    out_dir = os.path.dirname(junction_file_list)
    samples_tsv = os.path.join(out_dir, "samples.tsv")
    n_samples = 0

    with open(bam_list_path, "r") as f_in, \
         open(junction_file_list, "w") as f_junc_list, \
         open(samples_tsv, "w") as f_samples:

        reader = csv.reader(f_in, delimiter="\t")

        for row in reader:
            if not row or not row[0].strip() or row[0].lstrip().startswith("#"):
                continue

            bam_file = os.path.abspath(row[0].strip())
            sample_id = os.path.basename(bam_file).replace(".bam", "")
            group = row[1].strip() if len(row) > 1 and row[1].strip() else None
            donor = row[2].strip() if len(row) > 2 and row[2].strip() else None

            if not os.path.isfile(bam_file):
                logging.error(f"BAM not found, skipping: {bam_file}")
                continue

            junc_file = os.path.join(junction_output_dir, f"{sample_id}.junc")

            if os.path.exists(junc_file):
                logging.info(f"Junction file exists; skipping extraction: {junc_file}")
            else:
                if not os.path.exists(bam_file + ".bai"):
                    try:
                        pysam.index(bam_file)
                    except Exception as e:
                        logging.error(f"Failed to index BAM file {bam_file}: {e}")
                        continue

                logging.info(f"Extracting junctions from {bam_file} -> {junc_file}")
                cmd = [
                    "regtools", "junctions", "extract",
                    "-a", "6",
                    "-m", "10",
                    "-M", "1000000",
                    "-s", strand,
                    bam_file,
                    "-o", junc_file,
                ]
                subprocess.run(cmd, check=True)

            # Only reached when the .junc file exists, so the two lists stay in sync
            f_junc_list.write(junc_file + "\n")
            f_junc_list.flush()

            # Write samples.tsv: name, treatment, donor (donor = covariate in leafcutter_ds.R)
            if group:
                f_samples.write("\t".join(x for x in (sample_id, group, donor) if x) + "\n")
            else:
                f_samples.write(f"{sample_id}\n")
            f_samples.flush()
            n_samples += 1

    logging.info(f"Wrote: {junction_file_list}, {samples_tsv} ({n_samples} samples)")
    return junction_file_list, samples_tsv


def cluster_introns(junctions_file_list, out_dir, cell_type, script, py=sys.executable):
    os.makedirs(out_dir, exist_ok=True)

    # Leafcutter writes to CWD using prefix from -o
    cmd = [
        py, script,
        "-j", junctions_file_list,
        "-o", cell_type,
        "-m", "10",
        "-l", "1000000",
    ]
    subprocess.run(cmd, check=True, cwd=out_dir)

    # Tidy per-sample cluster files
    per_sample_clusters_dir = os.path.join(out_dir, "per_sample_clusters")
    os.makedirs(per_sample_clusters_dir, exist_ok=True)
    for file in glob.glob(os.path.join(out_dir, f"*.junc.{cell_type}.sorted.gz")):
        shutil.move(file, os.path.join(per_sample_clusters_dir, os.path.basename(file)))

    logging.info(f"Intron clustering finished. Prefix: {cell_type}")


def match_sample_names(counts_file, samples_tsv):
    """
    leafcutter_ds.R needs column 1 of samples.tsv to match the sample
    columns in the counts file exactly. Fix common suffix differences
    (e.g. a trailing .junc) and stop with a clear message otherwise.
    """
    with gzip.open(counts_file, "rt") as fh:
        header = fh.readline().split()
    count_names = set(header[1:] if header and header[0] == "chrom" else header)

    rows = [line.rstrip("\n").split("\t") for line in open(samples_tsv) if line.strip()]
    for r in rows:
        name = r[0]
        for candidate in (name, name + ".junc", name.replace(".junc", "")):
            if candidate in count_names:
                r[0] = candidate
                break
        else:
            logging.error(f"'{name}' is not a column in {counts_file}. "
                          f"Counts columns: {sorted(count_names)}")
            sys.exit(1)

    with open(samples_tsv, "w") as fh:
        for r in rows:
            fh.write("\t".join(r) + "\n")
    logging.info("Sample names in samples.tsv match the counts file")


def differential_quantifications(counts_file, group_file, out_dir, exon_file,
                                 min_samples_per_intron, cell_type, script_path,
                                 rscript="Rscript", threads=1):
    os.makedirs(out_dir, exist_ok=True)
    output_prefix = os.path.join(out_dir, cell_type)

    cmd = [
        rscript, script_path,
        counts_file,
        group_file,
        "-o", output_prefix,
        "-i", str(min_samples_per_intron),
        "-e", exon_file,
        "-p", str(threads),
    ]
    subprocess.run(cmd, check=True, cwd=out_dir)
    logging.info(f"Differential splicing finished. Results at prefix: {output_prefix}")


def main():
    parser = argparse.ArgumentParser(description="LeafCutter pipeline")
    parser.add_argument("--bam_file_list", type=str, required=True,
                        help="Headerless TSV: BAM path<TAB>treatment<TAB>donor")
    parser.add_argument("--output_dir", type=str, required=True, help="Output directory")
    parser.add_argument("--cell_type", type=str, required=True, help="Cell type label")
    parser.add_argument("--exon_file", type=str, required=True,
                        help="Exon .txt.gz from get_exon.sh / gtf_to_exons.R")
    parser.add_argument("--strand", default="RF", choices=["FR", "RF", "XS"],
                        help="regtools -s (default RF)")
    parser.add_argument("--leafcutter_dir", default=default_leafcutter_dir(),
                        help="Leafcutter clone (default: ./leafcutter or ../leafcutter next to this script)")
    parser.add_argument("--cluster_script", default=None,
                        help="Override path to leafcutter_cluster_regtools.py")
    parser.add_argument("--ds_script", default=None,
                        help="Override path to leafcutter_ds R script")
    parser.add_argument("--min_samples_per_intron", type=int, default=3)
    parser.add_argument("--threads", type=int,
                        default=int(os.environ.get("SLURM_CPUS_PER_TASK", "1")),
                        help="Threads for leafcutter_ds.R (default: SLURM_CPUS_PER_TASK or 1)")
    args = parser.parse_args()

    # ---- make every user path absolute NOW, relative to the launch folder,
    #      because the Leafcutter steps run with cwd=output_dir ----
    leafcutter_dir = os.path.abspath(args.leafcutter_dir)
    cluster_script = os.path.abspath(args.cluster_script or os.path.join(
        leafcutter_dir, "clustering", "leafcutter_cluster_regtools.py"))
    ds_script = os.path.abspath(args.ds_script or os.path.join(
        leafcutter_dir, "scripts", "leafcutter_ds.R"))
    bam_file_list = os.path.abspath(args.bam_file_list)
    exon_file = os.path.abspath(args.exon_file)
    output_dir = os.path.abspath(args.output_dir)

    # ---- logging FIRST ----
    os.makedirs(output_dir, exist_ok=True)
    log_dir = os.path.join(output_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, f"{args.cell_type}.log")
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s: %(message)s",
        handlers=[logging.FileHandler(log_file, mode="w"), logging.StreamHandler()],
        force=True,
    )
    logging.info(f"Logging to {log_file}")

    for f in (bam_file_list, exon_file, cluster_script, ds_script):
        if not os.path.isfile(f):
            logging.error(f"Not found: {f}")
            sys.exit(1)
    logging.info(f"Cluster script: {cluster_script}")
    logging.info(f"DS script:      {ds_script}")

    junction_output_dir = os.path.join(output_dir, "junction_files")
    junction_list_file = os.path.join(output_dir, "junction_files.txt")

    # ---- Step 1: Extract junctions + samples.tsv ----
    jlist, samples_tsv = extract_junctions(
        bam_file_list, junction_output_dir, junction_list_file, strand=args.strand
    )

    # ---- Step 2: Clustering ----
    cluster_introns(jlist, output_dir, args.cell_type, script=cluster_script)

    counts_file = os.path.join(output_dir, f"{args.cell_type}_perind_numers.counts.gz")
    if not os.path.isfile(counts_file):
        logging.error(f"Counts file not found: {counts_file}")
        raise FileNotFoundError(counts_file)

    match_sample_names(counts_file, samples_tsv)

    # ---- Step 3: Differential splicing (treatment, with donor as covariate) ----
    differential_quantifications(
        counts_file, samples_tsv, output_dir, exon_file,
        args.min_samples_per_intron, args.cell_type,
        script_path=ds_script, threads=args.threads
    )


if __name__ == "__main__":
    main()