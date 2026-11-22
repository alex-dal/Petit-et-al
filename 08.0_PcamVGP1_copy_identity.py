#!/usr/bin/env python3
"""Step 8a - pairwise nucleotide identity between copies of a feature.

Smith-Waterman (parasail, DNAfull), both orientations tested, identity taken over
the aligned columns. Sequences are upper-cased first: soft-masked assemblies
otherwise score as mismatches.

    python 08.0_PcamVGP1_copy_identity.py copies.fa copy_identity.json

copies.fa holds one record per locus, named by locus id (see 01.0/02.0).
"""
import sys,itertools,json,statistics,parasail,collections
def rd(p):
    d={};n=None;b=[]
    for l in open(p):
        if l[0]=='>':
            if n:d[n]=''.join(b).upper()
            n=l[1:].strip();b=[]
        else:b.append(l.strip())
    if n:d[n]=''.join(b).upper()
    return d
comp=str.maketrans('ACGTN','TGCAN'); rc=lambda s:s.translate(comp)[::-1]
seqs=rd(sys.argv[1]); names=sorted(seqs); mat=parasail.dnafull
def best(x,y):
    o=[]
    for z in (y,rc(y)):
        r=parasail.sw_trace_striped_16(x,z,10,1,mat); t=r.traceback
        al=len(t.comp)
        o.append((r.score, 100*t.comp.count('|')/al if al else 0, al))
    o.sort(reverse=True); return o[0]
res={}
for a,b in itertools.combinations(names,2):
    s,i,al=best(seqs[a],seqs[b]); res[(a,b)]=(i,al)
json.dump({f"{a}||{b}":v for (a,b),v in res.items()},open(sys.argv[2],'w'))
full=[n for n in names if len(seqs[n])>=5000]
sp=lambda n:n.split('_')[1]
def rep(v,lab):
    v=[x for x in v]
    if v: print("%-17s n=%3d median %.2f%%  range %.2f-%.2f"%(lab,len(v),statistics.median(v),min(v),max(v)))
rep([v[0] for (a,b),v in res.items() if a in full and b in full and sp(a)==sp(b)],"within-species")
rep([v[0] for (a,b),v in res.items() if a in full and b in full and sp(a)!=sp(b)],"between-species")
print("\nbetween-species identity distribution (full-length pairs):")
bins=collections.Counter()
for (a,b),v in res.items():
    if a in full and b in full and sp(a)!=sp(b): bins[int(v[0]//2)*2]+=1
for k in sorted(bins,reverse=True): print("   %3d-%3d%%  %s"%(k,k+2,'#'*bins[k]+' '+str(bins[k])))
