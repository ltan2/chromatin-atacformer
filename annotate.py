import snapatac2 as snap

def main():
    adata_file_names = ["GSM_D0.h5ad", "GSM_D2.h5ad", "GSM_D4.h5ad", "GSM_D6.h5ad", "GSM_D8.h5ad", "GSM_D10.h5ad", "GSM_D12.h5ad", "GSM_D14.h5ad"]

    for i, data_file in enumerate(adata_file_names):
        adata = snap.read(data_file), backed="r+")
        adata = adata.to_memory()   # <-- load everything into RAM

        snap.tl.macs3(
            adata,
            groupby=None,
            nolambda=True,
            shift=-100,
            extsize=200,
            key_added="peaks"
        )

        # Export peaks DataFrame to BED
        peak_bed_path = f"peaks_pseudobulk_{i+2}.bed"
        adata.uns["peaks_pseudobulk"].to_csv(
            peak_bed_path,
            sep="\t",
            header=False,
            index=False,
            columns=["chrom", "start", "end"]
        )

        # Use BED file path in make_peak_matrix
        snap.pp.make_peak_matrix(
            adata,
            peak_file=peak_bed_path,
            inplace=True
        )

        adata.write(f"peak_matrix_{i+2}.h5ad")

if __name__ == "__main__":
    main()
