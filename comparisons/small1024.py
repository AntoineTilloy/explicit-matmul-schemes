"""Draw the verified small-square improvements, keeping coefficient fields apart."""
import argparse
import csv
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--field',choices=('Q','F2','all'),default='all')
    args=parser.parse_args()
    rational=json.loads((ROOT/'checks/small1024_rectangular.json').read_text())['changes']
    rational=[dict(field='Q',n=p['n'],rank=p['rank'],baseline_rank=p['comparison_rank'],
                   construction=p['classification']) for p in rational]
    asymmetric=json.loads((ROOT/'checks/asymmetric_strassen.json').read_text())['changes']
    rational.extend(dict(field='Q',n=p['n'],rank=p['rank'],baseline_rank=p['comparison_rank'],
                         construction=p['classification']) for p in asymmetric)
    gauge=json.loads((ROOT/'checks/small1024_recompression.json').read_text())['changes'][0]
    rational.append(dict(field='Q',n=1024,rank=int(gauge['rank']),
                         baseline_rank=int(gauge['baseline_rank']),construction='rational CW recompression, odd real character and square cancellation'))
    binary=[]
    for n in (61,243,969):
        root=ROOT/'characteristic2'
        spec=json.loads((root/f'F2_{n}_class_recompressed.json').read_text())
        base=json.loads((root/spec['base']).read_text())
        binary.append(dict(field='F2',n=n,rank=spec['rank'],baseline_rank=base['baseline'],
                           construction='permuted restriction and exact pair merges'))
    points=sorted(rational,key=lambda p:p['n'])+binary
    for p in points:
        p['saving']=p['baseline_rank']-p['rank'];assert p['saving']>0
        p['saving_percent']=100*p['saving']/p['baseline_rank']
    with (ROOT/'comparisons/small1024.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=points[0].keys());writer.writeheader();writer.writerows(points)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none',
                         'axes.spines.top':False,'axes.spines.right':False})
    for field,color,title,destination in (
            ('Q','#176c99','Square improvements over R and C',ROOT/'figures'/'small1024'),
            ('F2','#348568','Square improvements over F2',ROOT/'characteristic2'/'improvements')):
        if args.field!='all' and args.field!=field:continue
        fig,ax=plt.subplots(figsize=(9.2,8.2 if field=='Q' else 5.8))
        rows=[p for p in points if p['field']==field];x=list(range(len(rows)))
        ax.barh(x,[p['saving_percent'] for p in rows],color=color,height=.60,zorder=3)
        for j,p in enumerate(rows):
            percent=f"{p['saving_percent']:.6f}%" if p['n']==1024 else f"{p['saving_percent']:.2f}%"
            ax.text(p['saving_percent']+.035,j,percent+' · '+f"{p['rank']:,} products",
                    ha='left',va='center',fontsize=8.5,color=color)
        ax.set_yticks(x,[str(p['n'])+(' (LITA refresh)' if p['n']==432 else '') for p in rows])
        ax.set_ylabel('Square matrix size n',labelpad=9)
        ax.set_ylim(len(rows)-.4,-.6);ax.set_xlim(0,3.0);ax.grid(axis='x',color='#edf0f3');ax.set_axisbelow(True)
        ax.set_xlabel('Products saved against the stated reference (%)')
        fig.suptitle(title,x=.10,y=.98,ha='left',fontsize=18,weight='bold')
        subtitle=('Exact rational coefficients • Square sizes ≤ 1,024' if field=='Q'
                  else 'Exact F2 coefficients • Separate from the R/C catalogue')
        fig.text(.10,.91,subtitle,color='#526171',fontsize=10.5)
        footer=('References: previous strengthened composition envelope / previous saved 432 and 1024.\n'
                'The 432 change refreshes LITA primitives. Sub-2.7 at size 1,024 remains open.' if field=='Q'
                else 'Reference: strengthened finite F2 composition closure.\n'
                'These schemes do not establish corresponding counts over R or C.')
        fig.text(.10,.035,footer,fontsize=9,color='#526171')
        fig.subplots_adjust(left=.22,right=.985,top=.83,bottom=.15)
        for suffix in ('svg','png','pdf'):
            fig.savefig(destination.with_suffix('.'+suffix),dpi=170,facecolor='white')
        plt.close(fig)
    print(json.dumps(dict(points=points,scope='Exact compositional constructions with stated primitive proof basis')))


if __name__=='__main__':main()
