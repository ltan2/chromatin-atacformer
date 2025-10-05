import scanpy as sc
from geniml.atacformer import AtacformerForCellClustering
from gtars.tokenizers import Tokenizer, tokenize_fragment_file
from geniml.tokenization.utils import tokenize_anndata
import numpy as np
import glob

def main():
    atac_model = AtacformerForCellClustering.from_pretrained("databio/atacformer-base-hg38")
    tokenizer = Tokenizer.from_pretrained("databio/atacformer-base-hg38")
    

    files = glob.glob("peak_matrix_*.h5ad")

    # read them into a list of AnnData objects
    adatas = [sc.read(f) for f in files]
    
    for i,data in enumerate(adatas):

        peaks = data.var_names.to_series().str.split("[:-]", expand=True)
        peaks.columns = ["chr", "start", "end"]

        data.var = peaks
        data.var["start"] = data.var["start"].astype(int)
        data.var["end"] = data.var["end"].astype(int)

        batch_size = 32
        embeddings_list = []
        for i in range(0, len(data), batch_size):
            batch = data[i:i+batch_size]
            tokens = tokenize_anndata(batch, tokenizer)

            input_ids = [t["input_ids"] for t in tokens]
            embeddings = atac_model.encode_tokenized_cells(
                input_ids=input_ids,
                batch_size=2
            )
            
            embeddings_list.append(embeddings.detach().cpu().numpy())

        all_embeddings = np.vstack(embeddings_list)

        # Save to disk
        np.save(f"embeddings_d{i+2}.npy", all_embeddings)  # .npy format

if __name__ == "__main__":
    main()
