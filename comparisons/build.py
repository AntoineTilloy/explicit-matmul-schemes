"""Recompute finite bounds and plots from pinned primitives and saved recipes."""
import csv
from fractions import Fraction as F
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent


def lita(n):
    return int(F(n**3,3)+3*n*n+F(37*n,6)+5)


def close(limit,seeds,extra=()):
    ranks=[n**3 for n in range(limit+1)]
    recipes=[dict(kind='classical',n=n) for n in range(limit+1)]
    for n in range(1,limit+1):
        options=[]
        if str(n) in seeds:options.append((int(seeds[str(n)]['rank']),dict(kind='published',n=n,source=str(n))))
        if n>=8 and n%2==0:options.append((lita(n),dict(kind='lita-even',n=n)))
        if 9<=n<32 and n%2:options.append((int(F(n**3,3)+F(7*n*n,2)+F(14*n,3)-F(9,2)),dict(kind='lita-odd',n=n)))
        for r,why in options:
            if r<ranks[n]:ranks[n],recipes[n]=r,why
    for entry in extra:
        n=int(entry['n']);r=int(entry['rank'])
        if n<=limit and r<ranks[n]:ranks[n],recipes[n]=r,dict(kind='catalogue',n=n,rank=str(r),root_object=entry['root_object'])
    for iteration in range(100):
        changed=False
        for n in range(2,limit+1):
            r=ranks[n-1]+3*n*n-3*n+1
            if r<ranks[n]:ranks[n],recipes[n]=r,dict(kind='peel',n=n,child=recipes[n-1]);changed=True
            for a in range(2,math.isqrt(n)+1):
                if n%a==0 and ranks[a]*ranks[n//a]<ranks[n]:
                    ranks[n]=ranks[a]*ranks[n//a];recipes[n]=dict(kind='product',n=n,left=recipes[a],right=recipes[n//a]);changed=True
        for n in range(limit-1,0,-1):
            if ranks[n+1]<ranks[n]:ranks[n]=ranks[n+1];recipes[n]=dict(kind='pad',n=n,child=recipes[n+1]);changed=True
        if not changed:return ranks,recipes
    raise RuntimeError('Closure did not converge')


def evaluate(d,seeds):
    k=d['kind'];n=d['n']
    if k=='classical':return n**3
    if k in ('lita-even','LITA'):return lita(n)
    if k=='lita-odd':return int(F(n**3,3)+F(7*n*n,2)+F(14*n,3)-F(9,2))
    if k=='published':return int(seeds[d['source']]['rank'])
    if k in ('saved-comparison-primitive','catalogue'):return int(d['rank'])
    if k=='product':
        assert n==d['left']['n']*d['right']['n'];return evaluate(d['left'],seeds)*evaluate(d['right'],seeds)
    if k=='pad':assert n<=d['child']['n'];return evaluate(d['child'],seeds)
    if k=='peel':assert n>d['child']['n'];return evaluate(d['child'],seeds)+n**3-d['child']['n']**3
    raise ValueError(k)


def main():
    entries=json.loads((ROOT/'catalogue.json').read_text());seeds=json.loads((HERE/'primitives.json').read_text())
    ranks,recipes=close(16384,seeds);aug,_=close(16384,seeds,entries)
    large=json.loads((HERE/'large_recipes.json').read_text());points=[];witnesses={}
    for entry in entries:
        n=int(entry['n']);rank=int(entry['rank'])
        recipe=recipes[n] if n<=16384 else large[str(n)]
        comparison=evaluate(recipe,seeds)
        assert n>16384 or comparison==ranks[n]
        witnesses[str(n)]=recipe
        points.append(dict(n=n,rank=rank,effective_exponent=math.log(rank)/math.log(n),
                           comparison_rank=comparison,comparison_exponent=math.log(comparison)/math.log(n),
                           scope='finite closure' if n<=16384 else 'specific published construction',
                           improves_comparison=rank<comparison))
    with (HERE/'counts.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(points[0]));w.writeheader();w.writerows(points)
    (HERE/'recipes.json').write_text(json.dumps(witnesses,separators=(',',':'))+'\n')
    with (HERE/'closure.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['n','literature_rank','with_catalogue_rank'])
        w.writerows((n,ranks[n],aug[n]) for n in range(2,16385))
    draw(points,ranks,aug)
    print(json.dumps({'points':len(points),'losing_sizes':[p['n'] for p in points if not p['improves_comparison']],'lita32_rank':lita(32)}))


def draw(points,ranks,aug):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(14,6.3),gridspec_kw={'width_ratios':[1.1,1]})
    blue='#176c99';orange='#c17529';ref=math.log(lita(32))/math.log(32)
    ns=list(range(400,16385))
    axes[0].plot(ns,[math.log(ranks[n])/math.log(n) for n in ns],color=orange,lw=.75,alpha=.65,label='Pinned literature + finite closure')
    axes[0].plot(ns,[math.log(aug[n])/math.log(n) for n in ns],color=blue,lw=.85,alpha=.8,label='Literature + catalogue + same closure')
    for ax,subset,title in [(axes[0],[p for p in points if p['n']<=16384],'Detail through 16,384'),(axes[1],points,'All saved sizes')]:
        for p in subset:
            ax.plot([p['n']]*2,[p['effective_exponent'],p['comparison_exponent']],color='#d6dce2',lw=.8,zorder=1)
        win=[p for p in subset if p['improves_comparison']];loss=[p for p in subset if not p['improves_comparison']]
        ax.scatter([p['n'] for p in win],[p['effective_exponent'] for p in win],color=blue,s=25,zorder=4,label='Saved scheme')
        ax.scatter([p['n'] for p in loss],[p['effective_exponent'] for p in loss],color='#b74349',marker='x',s=32,zorder=5,label='Saved candidate loses comparison')
        ax.axhline(ref,color='#8561a2',ls='--',lw=1,label=f'LITA32 reference: {ref:.6f}')
        ax.set_xscale('log',base=2);ax.set_title(title,loc='left',pad=12);ax.grid(axis='y',color='#edf0f3');ax.set_axisbelow(True)
        ax.set_xlabel('Square matrix size n (log scale)',labelpad=10)
    large=[p for p in points if p['n']>16384];small=[p for p in points if p['n']<=16384]
    axes[1].scatter([p['n'] for p in small],[p['comparison_exponent'] for p in small],marker='s',s=25,facecolors='white',edgecolors=orange,label='Finite literature closure')
    axes[1].scatter([p['n'] for p in large],[p['comparison_exponent'] for p in large],marker='D',s=28,color=orange,label='Specific published construction')
    axes[0].set_ylabel(r'Effective exponent $\log_n R(n)$ — lower is better')
    axes[0].set_xlim(390,18000);axes[0].set_ylim(2.726,2.79)
    axes[0].set_xticks([512,1024,2048,4096,8192,16384]);axes[0].xaxis.set_major_formatter(FuncFormatter(lambda x,pos:f'{int(x):,}'))
    powers=[9,14,19,24,29,34];axes[1].set_xticks([2**k for k in powers],[f'$2^{{{k}}}$' for k in powers]);axes[1].set_xlim(330,3e10);axes[1].set_ylim(2.69,2.79)
    axes[1].axhline(2.7,color='#479b7d',ls=':',lw=.9)
    axes[1].annotate('2.698590 at n ≈ 13.47 billion',(13468840704,2.6985901148),xytext=(-180,-15),textcoords='offset points',fontsize=9,arrowprops={'arrowstyle':'-','color':'#a4adb5'})
    fig.suptitle('Explicit square matrix multiplication schemes',x=.075,y=.98,ha='left',fontsize=19,weight='bold')
    fig.text(.075,.916,'Exact bilinear counts over Q • Best available bounds in this catalogue • Pinned literature comparison',color='#526171',fontsize=10.5)
    h,l=axes[0].get_legend_handles_labels();h2,l2=axes[1].get_legend_handles_labels()
    for hh,ll in zip(h2,l2):
        if ll not in l:h.append(hh);l.append(ll)
    fig.legend(h,l,loc='lower left',bbox_to_anchor=(.07,.045),ncol=3,frameon=False,fontsize=8.5)
    fig.text(.075,.018,'Dense curves use the same declared composition closure; large diamonds are individual constructions. No claim of globally optimal ranks.',fontsize=8.5,color='#526171')
    fig.subplots_adjust(left=.075,right=.98,top=.825,bottom=.265,wspace=.22)
    (ROOT/'figures').mkdir(exist_ok=True)
    for suffix in ('svg','png','pdf'):fig.savefig(ROOT/'figures'/f'exponents.{suffix}',dpi=170,facecolor='white')


if __name__=='__main__':main()
