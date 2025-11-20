import pandas as pd
import mygene

# ----------------------------
# INPUT: Your gene list
# ----------------------------
# Example: load the gene panel you extracted
df = pd.read_csv("datasets_mouse_brain_map_BrainReceptorShowcase_Slice1_Replicate1_cell_by_gene_S1R1.csv")
genes = [c for c in df.columns[1:] if not c.lower().startswith("unnamed")]

genes = list(dict.fromkeys(genes))  # remove duplicates, preserve order

# ----------------------------
# QUERY ENSEMBL VIA mygene
# ----------------------------
mg = mygene.MyGeneInfo()

# Query Ensembl transcript IDs for mouse
out = mg.querymany(
    genes,
    scopes='symbol',
    fields='ensembl.transcript',
    species='mouse'
)

# ----------------------------
# FORMAT RESULTS
# ----------------------------
rows = []
for q in out:
    symbol = q['query']
    if 'ensembl' in q:
        ens = q['ensembl']
        # Sometimes ensembl field is a list
        if isinstance(ens, list):
            transcripts = [e.get('transcript') for e in ens if 'transcript' in e]
            transcript_id = transcripts[0] if transcripts else None
        else:
            transcript_id = ens.get('transcript')
    else:
        transcript_id = None
    
    rows.append({
        'gene_symbol': symbol,
        'transcript_id': transcript_id
    })

mapped = pd.DataFrame(rows)

# ----------------------------
# SAVE TO FILE
# ----------------------------
mapped.to_csv("mouse_panel_transcript_ids.csv", index=False)
print("Saved mouse_panel_transcript_ids.csv with", len(mapped), "genes.")
