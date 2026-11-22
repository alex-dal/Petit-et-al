#!/usr/bin/env python3
"""Step 8b - extract the protein-coding genes flanking each locus.

Genes overlapping a 50 kb window on either side of a locus, and not overlapping
the element itself, are taken from the annotation and translated using the CDS
phase. Peptides shorter than 60 aa or carrying an internal stop are dropped.
A window is marked complete only when 50 kb fits inside the contig on both sides:
a window truncated by a contig end yields fewer genes and an unstable denominator.

Augustus annotations may name contigs scaffold_N; scaffold_N is taken to be the
Nth longest contig of that assembly.

    python 08.1_PcamVGP1_flanking_genes.py loci.json annotations/ genomes/ flanking.json
"""
import sys,os,json,re,collections
LOCI,GFFDIR,FADIR,OUT=sys.argv[1:5]
W=50000
loci=[L for L in json.load(open(LOCI)) if 'id' in L]
gffs={}
for f in os.listdir(GFFDIR):
    m=re.match(r'(GC[AF]_\d+\.\d+)',f)
    if m and (f.endswith('.gff') or f.endswith('.gff3')):
        gffs.setdefault(m.group(1),[]).append(f)
def pick(acc):
    fs=gffs.get(acc,[])
    for f in fs:
        if f.endswith('.gff'): return f
    return fs[0] if fs else None
def contig_order(acc):
    d=os.path.join(FADIR,acc); fn=[x for x in os.listdir(d) if x.endswith('.fna')][0]
    lens={};name=None;n=0
    for l in open(os.path.join(d,fn)):
        if l[0]=='>':
            if name:lens[name]=n
            name=l[1:].split()[0];n=0
        else:n+=len(l.strip())
    if name:lens[name]=n
    return lens,[c for c,_ in sorted(lens.items(),key=lambda x:-x[1])]
def load_cds(acc):
    f=pick(acc)
    if not f: return None,None
    lens,order=contig_order(acc)
    scaf = 'scaffold_' in open(os.path.join(GFFDIR,f)).readline(4000) or True
    cds=collections.defaultdict(list); gmap={}
    for line in open(os.path.join(GFFDIR,f)):
        if line[0]=='#': continue
        p=line.rstrip('\n').split('\t')
        if len(p)<9 or p[2]!='CDS': continue
        sid=p[0]
        if sid.startswith('scaffold_'):
            i=int(sid.split('_')[1])-1
            if i>=len(order): continue
            sid=order[i]
        par=re.search(r'Parent=([^;]+)',p[8])
        key=(sid,par.group(1) if par else p[8][:30])
        cds[key].append((int(p[3]),int(p[4]),p[6],p[7]))
    return cds,lens
comp=str.maketrans('ACGTNacgtn','TGCANtgcan')
CODON={}
B='TCAG'
AAs='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
i=0
for a in B:
    for b in B:
        for c in B: CODON[a+b+c]=AAs[i]; i+=1
def translate(s):
    s=s.upper(); return ''.join(CODON.get(s[i:i+3],'X') for i in range(0,len(s)-2,3))
def getseq(acc,contig):
    d=os.path.join(FADIR,acc); fn=[x for x in os.listdir(d) if x.endswith('.fna')][0]
    cur=None;buf=[]
    for l in open(os.path.join(d,fn)):
        if l[0]=='>':
            if cur==contig: return ''.join(buf)
            cur=l[1:].split()[0];buf=[]
        elif cur==contig: buf.append(l.strip())
    return ''.join(buf) if cur==contig else None
cache={}
out={}
for L in loci:
    acc=L['acc']
    if acc not in cache: cache[acc]=load_cds(acc)
    cds,lens=cache[acc]
    if cds is None: out[L['id']]=None; continue
    lo,hi=L['start']-W,L['end']+W
    genes=[]
    for (sid,par),parts in cds.items():
        if sid!=L['contig']: continue
        a=min(x[0] for x in parts); b=max(x[1] for x in parts)
        if b<lo or a>hi: continue
        if not (b < L['start'] or a > L['end']):  # overlaps element itself
            continue
        genes.append((a,b,par,parts))
    seq=getseq(acc,L['contig'])
    prots={}; pos={}
    for a,b,par,parts in genes:
        parts=sorted(parts)
        s=''.join(seq[x[0]-1:x[1]] for x in parts)
        if parts[0][2]=='-': s=s.translate(comp)[::-1]
        try: ph=int(parts[0][3])
        except: ph=0
        p=translate(s[ph:]).rstrip('*')
        if len(p)>=60 and '*' not in p[:-1]: prots[par]=p; pos.setdefault(par,(a,b,parts[0][2]))
    out[L['id']]=dict(n_genes=len(genes),n_prot=len(prots),prots=prots,pos={k:list(v) for k,v in pos.items()},
                      contig_len=L.get('contig_len'),dist_end=L.get('dist_to_end'),
                      species=L['species'],strain=L['strain'],acc=acc,
                      window_full=(L['start']-W>=0 and L['end']+W<=L.get('contig_len',0)))
json.dump(out,open(OUT,'w'))
ok=[k for k,v in out.items() if v]
print("loci with annotation: %d / %d"%(len(ok),len(loci)))
print("loci with full +/-50kb window: %d"%sum(1 for k in ok if out[k]['window_full']))
import statistics
ns=[out[k]['n_prot'] for k in ok]
print("flanking proteins per locus: median %d  range %d-%d  total %d"%(statistics.median(ns),min(ns),max(ns),sum(ns)))
for k in sorted(ok):
    v=out[k]
    print("  %-50s genes=%3d prot=%3d ctg=%9d distend=%7d %s"%(k[:50],v['n_genes'],v['n_prot'],v['contig_len'] or 0,v['dist_end'] or 0,"" if v['window_full'] else "TRUNCATED"))
