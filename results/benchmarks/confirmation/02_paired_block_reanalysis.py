"""Replicate-block bootstrap, preserving cross-cell and method dependence.

Inputs are existing per-replication summaries, not individual coefficients.
Complete replicate blocks are required for all cells in each comparison.
"""
import collections
import csv
import hashlib
import math
import os
import pathlib
import random
import statistics

# Defaults are relative to the repository root. The paper projects keep the raw confirmation CSVs in their own
# results/benchmarks_confirmation/ folder, so the input and output folders can be overridden by environment variable.
ROOT = pathlib.Path(os.environ.get('SYDNI_CONFIRMATION_DIR', 'results/benchmarks/confirmation'))
OUT = pathlib.Path(os.environ.get('SYDNI_AUDIT_OUT', 'validation/audits/02_reanalysis'))
OUT.mkdir(parents=True, exist_ok=True)

def read(name):
    p = ROOT / name
    with p.open(newline='', encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))

def quantile(values, p):
    s=sorted(values); x=(len(s)-1)*p; lo=math.floor(x); hi=math.ceil(x)
    return s[lo]+(x-lo)*(s[hi]-s[lo])

def block_pair(rows,left,right,cells,metric,reps,seed=20261105,boot=5000):
    index={}
    for row in rows:
        if row['method'] not in (left,right): continue
        cell=(int(row['n_obs']),float(row['aug_ratio']))
        if cell not in cells: continue
        key=(int(row['replicate']),cell,row['method'])
        if key in index: raise ValueError('duplicate key: '+str(key))
        index[key]=row
    blocks=[];incomplete=[];cell_counts=collections.Counter()
    for rep in range(1,reps+1):
        differences=[]
        for cell in cells:
            a=index.get((rep,cell,left));b=index.get((rep,cell,right))
            if a is None or b is None or a['failed']!='FALSE' or b['failed']!='FALSE': continue
            try: x=float(a[metric]);y=float(b[metric])
            except (ValueError,KeyError): continue
            if not math.isfinite(x) or not math.isfinite(y): continue
            cell_counts[cell]+=1;differences.append(x-y)
        if len(differences)==len(cells): blocks.append(statistics.mean(differences))
        else: incomplete.append(rep)
    if len(blocks)<2: raise ValueError('too few complete replicate blocks')
    rng=random.Random(seed); draws=[]
    for _ in range(boot): draws.append(statistics.mean(rng.choices(blocks,k=len(blocks))))
    return dict(left=left,right=right,metric=metric,planned_blocks=reps,
        complete_blocks=len(blocks),excluded_blocks=len(incomplete),
        excluded_ids=';'.join(map(str,incomplete)),matched_pairs=sum(cell_counts.values()),
        min_matched_per_cell=min(cell_counts.values()),max_matched_per_cell=max(cell_counts.values()),
        estimate=statistics.mean(blocks),lower=quantile(draws,.025),upper=quantile(draws,.975),
        bootstrap_replicates=boot,seed=seed)

# A constant difference must be preserved exactly by the bootstrap.
fixture=[]
for rep in range(1,4):
 for n in (75,150):
  for method,value in (('a',3),('b',1)):
   fixture.append(dict(replicate=str(rep),n_obs=str(n),aug_ratio='1',method=method,failed='FALSE',metric=str(value)))
test=block_pair(fixture,'a','b',[(75,1.),(150,1.)],'metric',3,boot=20)
assert test['estimate']==test['lower']==test['upper']==2

s1=[r for r in read('study1_nonbrm_confirmation_raw.csv') if r['scenario']=='null_only']
cells1=[(n,ratio) for n in (75,150,300) for ratio in (.5,1.,2.)]
results=[]
for method in ('sydni_copula','sydni_graphical','sydni_bart'):
 for baseline in ('observed_only','bootstrap'):
  z=block_pair(s1,method,baseline,cells1,'type1_rate',100)
  z['study']='study1_null_only';z['grid']='9 cells, n=75/150/300, ratio=.5/1/2';results.append(z)

s2=read('study2_nonbrm_confirmation_raw.csv')+read('study2_brm_confirmation_raw.csv')
cells2=[(n,ratio) for n in (75,150) for ratio in (1.,2.)]
for engine in ('copula','graphical','bart','factor','brm'):
 for metric in ('within_scale_cor_error','within_scale_alpha_abs_change','between_construct_magnitude_change'):
  z=block_pair(s2,'sydni_'+engine+'_ma','sydni_'+engine+'_flat',cells2,metric,50 if engine=='brm' else 100)
  z['study']='study2_ma_minus_flat';z['grid']='4 cells, n=75/150, ratio=1/2';results.append(z)
fields=['study','grid','left','right','metric','planned_blocks','complete_blocks','excluded_blocks','excluded_ids','matched_pairs','min_matched_per_cell','max_matched_per_cell','estimate','lower','upper','bootstrap_replicates','seed']
with (OUT/'paired_block_bootstrap.csv').open('w',newline='',encoding='utf-8') as f:
 writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(results)
for z in results:
 print(z['study'],z['left'],'minus',z['right'],z['metric'],'blocks',z['complete_blocks'],'/',z['planned_blocks'],
       'estimate',round(z['estimate'],5),'interval',round(z['lower'],5),round(z['upper'],5))
for name in ('study1_nonbrm_confirmation_raw.csv','study2_nonbrm_confirmation_raw.csv','study2_brm_confirmation_raw.csv'):
 print('sha256',name,hashlib.sha256((ROOT/name).read_bytes()).hexdigest())
