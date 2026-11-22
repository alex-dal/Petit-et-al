#!/usr/bin/env python3
"""Step 8c - orthologues between the flanking windows of two loci.

Proteins are compared by local alignment (parasail, BLOSUM62) after a 4-mer
prefilter, and orthologues are reciprocal best hits covering at least half the
longer protein. Conservation between two loci is the number of orthologous pairs
divided by the smaller of the two gene counts, reported only where both windows
are complete.

    python 08.2_PcamVGP1_flanking_orthologs.py flanking.json copy_identity.json conservation.json
"""
import sys,json,itertools,statistics,collections,parasail
FL,IDENT,OUT=sys.argv[1:4]
fl=json.load(open(FL))
loci=[k for k in sorted(fl) if fl[k] and fl[k]['n_prot']>0]
P=[];owner=[]
for k in loci:
    for name,seq in fl[k]['prots'].items():
        P.append(seq.upper()); owner.append(k)
print("proteins:",len(P),"across",len(loci),"loci")
mat=parasail.blosum62
def kset(s,k=4): return {s[i:i+k] for i in range(len(s)-k+1)}
K=[kset(s) for s in P]
best=collections.defaultdict(lambda:(-1,None))   # (score, j) best hit of i within locus of j
pairsc={}
n=len(P)
for i in range(n):
    for j in range(i+1,n):
        if owner[i]==owner[j]: continue
        a,b=K[i],K[j]
        if len(a&b)/min(len(a),len(b)) < 0.10: continue
        r=parasail.sw_striped_16(P[i],P[j],11,1,mat)
        sc=r.score
        cov=min(len(P[i]),len(P[j]))/max(len(P[i]),len(P[j]))
        if cov<0.5: continue
        norm=sc/ (3.0*min(len(P[i]),len(P[j])))
        if norm<0.35: continue
        pairsc[(i,j)]=sc
        if sc>best[(i,owner[j])][0]: best[(i,owner[j])]=(sc,j)
        if sc>best[(j,owner[i])][0]: best[(j,owner[i])]=(sc,i)
rbh=set()
for (i,lj),(sc,j) in best.items():
    if j is None: continue
    sc2,i2=best.get((j,owner[i]),(-1,None))
    if i2==i: rbh.add((min(i,j),max(i,j)))
print("reciprocal best-hit ortholog pairs:",len(rbh))
cnt=collections.Counter()
for i,j in rbh: cnt[(owner[i],owner[j])]+=1; cnt[(owner[j],owner[i])]+=1
res={}
for a,b in itertools.combinations(loci,2):
    na,nb=fl[a]['n_prot'],fl[b]['n_prot']
    sh=cnt.get((a,b),0)
    res[f"{a}||{b}"]=dict(shared=sh,frac=sh/min(na,nb) if min(na,nb) else 0,na=na,nb=nb,
                          full=fl[a]['window_full'] and fl[b]['window_full'])
json.dump(res,open(OUT,'w'))
sp=lambda k: fl[k]['species']
def rep(sel,lab):
    v=[x for x in sel]
    if v: print("%-34s n=%3d median %.1f%%  mean %.1f%%"%(lab,len(v),100*statistics.median(v),100*statistics.mean(v)))
full=lambda a,b: fl[a]['window_full'] and fl[b]['window_full']
rep([v['frac'] for k,v in res.items() for a,b in [k.split('||')] if sp(a)==sp(b) and full(a,b)],"within-species (full windows)")
rep([v['frac'] for k,v in res.items() for a,b in [k.split('||')] if sp(a)!=sp(b) and full(a,b)],"between-species (full windows)")
rep([v['frac'] for k,v in res.items() for a,b in [k.split('||')] if sp(a)!=sp(b) and not full(a,b)],"between-species (truncated)")
print("\nhighest flanking conservation between species (full windows):")
tops=sorted(((v['frac'],k) for k,v in res.items() for a,b in [k.split('||')] if sp(a)!=sp(b) and full(a,b)),reverse=True)[:8]
for f,k in tops:
    a,b=k.split('||'); print("   %5.1f%%  %s  vs  %s"%(100*f,a[:34],b[:34]))
