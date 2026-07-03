import pandas as pd
import os
import glob

def merge_overlaps(df):
    merged_rows = []

    # Group by Query, strand, frame, Accession
    for (query, strand, frame, acc), group in df.groupby(['Query', 'strand', 'frame', 'Accession']):
        # Sort by genomic start
        group = group.sort_values('start').reset_index(drop=True)

        current_start = None
        current_end = None
        current_rows = []

        for idx, row in group.iterrows():
            row_start = row['start']
            row_end = row['end']

            if current_start is None:
                current_start = row_start
                current_end = row_end
                current_rows = [row]
            elif row_start <= current_end:  # Overlaps
                current_end = max(current_end, row_end)
                current_rows.append(row)
            else:
                # Merge previous block using the best E-value
                best_hit = min(current_rows, key=lambda r: r['E-Value'])
                merged_row = best_hit.copy()
                merged_row['start'] = current_start
                merged_row['end'] = current_end

                # Compute DNA coordinates of protein
                merged_row['startp'] = merged_row['start'] + (best_hit['From'] - 1) * 3
                merged_row['endp']   = merged_row['start'] + (best_hit['To'] - 1) * 3

                merged_rows.append(merged_row)

                # Start new block
                current_start = row_start
                current_end = row_end
                current_rows = [row]

        # Handle last block
        if current_rows:
            best_hit = min(current_rows, key=lambda r: r['E-Value'])
            merged_row = best_hit.copy()
            merged_row['start'] = current_start
            merged_row['end'] = current_end
            merged_row['startp'] = merged_row['start'] + (best_hit['From'] - 1) * 3
            merged_row['endp']   = merged_row['start'] + (best_hit['To'] - 1) * 3
            merged_rows.append(merged_row)

    return pd.DataFrame(merged_rows)

def remove_overlapping_domains(df):
    # Sort by startp
    df = df.sort_values(['Query', 'startp', 'endp']).reset_index(drop=True)
    to_keep = []

    for query, group in df.groupby('Query'):
        group = group.sort_values('startp').reset_index(drop=True)
        prev = None
        for idx, row in group.iterrows():
            if prev is None:
                prev = row
                to_keep.append(row)
            else:
                # Check overlap
                if row['startp'] <= prev['endp']:  # overlapping
                    len_prev = prev['To'] - prev['From']
                    len_row = row['To'] - row['From']
                    if len_row > len_prev:
                        # Replace prev with row
                        to_keep[-1] = row
                        prev = row
                    else:
                        # Keep prev, skip row
                        continue
                else:
                    prev = row
                    to_keep.append(row)

    return pd.DataFrame(to_keep)

def process_file(filepath):
    df = pd.read_csv(filepath, sep="\t")
    merged_df = merge_overlaps(df)
    filtered_df = remove_overlapping_domains(merged_df)

    outpath = os.path.basename(filepath).replace("_cln.txt", "_merged_filtered.txt")
    filtered_df.to_csv(outpath, sep="\t", index=False)
    print(f"Merged and filtered file saved to {outpath}")

if __name__ == "__main__":
    input_files = glob.glob("*_cln.txt")
    for file in input_files:
        process_file(file)