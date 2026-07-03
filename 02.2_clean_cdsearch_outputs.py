import pandas as pd
import re
import os
import sys

def clean_query(x):
    # Remove leading Q#N - 
    x = re.sub(r"^Q#\d+\s*-\s*", "", x)
    # Keep only the part after the last ::
    x = x.split("::")[-1]
    # Remove the _strand*_frame*_start_end part at the end
    x = re.sub(r"_strand[-+]?\d+_frame\d+_\d+_\d+$", "", x)
    return x

def parse_last_field(x):
    # Example: strand-1_frame2_67427_70816
    m = re.search(r"_strand([-+]?\d+)_frame(\d+)_([0-9]+)_([0-9]+)$", x)
    if m:
        return pd.Series({
            'strand': m.group(1),
            'frame': m.group(2),
            'start': m.group(3),
            'end': m.group(4)
        })
    else:
        return pd.Series({'strand': None, 'frame': None, 'start': None, 'end': None})

def simplify_file(filepath):
    df = pd.read_csv(filepath, sep="\t")
    
    # Save the original first column for parsing
    orig_col = df.iloc[:,0]
    
    # Parse strand/frame/start/end
    df[['strand','frame','start','end']] = orig_col.apply(parse_last_field)
    
    # Clean first column
    df.iloc[:,0] = orig_col.apply(clean_query)
    
    # Save simplified output
    dirpath = os.path.dirname(filepath)
    base = os.path.splitext(os.path.basename(filepath))[0]
    outpath = os.path.join(dirpath, f"{base}_cln.txt")
    df.to_csv(outpath, sep="\t", index=False)
    print(f"Saved simplified file to {outpath}")

if __name__ == "__main__":
    for file in sys.argv[1:]:
        if os.path.isfile(file):
            simplify_file(file)
        else:
            print(f"Warning: {file} is not a valid file, skipping.")