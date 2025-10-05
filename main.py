from geniml.atacformer import AtacformerForCellClustering
import snapatac2 as snap

def main():
    fragment_files = [
        "GSM7763395_D0.frag.bed.gz",
        "GSM7763396_D2.frag.bed.gz",
        "GSM7763397_D4.frag.bed.gz",
        "GSM7763398_D6.frag.bed.gz",
        "GSM7763399_D8.frag.bed.gz",
        "GSM7763400_D10.frag.bed.gz",
        "GSM7763401_D12.frag.bed.gz",
        "GSM7763402_D14.frag.bed.gz"
    ]

    data_objects = [
        snap.pp.import_fragments(file, chrom_sizes=snap.genome.hg38, sorted_by_barcode=False)
        for file in fragment_files
    ]

    for i,data in enumerate(data_objects):
        data = data_objects[i]
        snap.metrics.tsse(data, snap.genome.hg38)
        snap.pl.tsse(data, interactive=True)
        snap.pp.filter_cells(data, min_counts=12000, min_tsse=15, max_counts=30000)

        data.write(f"GSM_D{i+2}.h5ad")


if __name__ == "__main__":
    main()



