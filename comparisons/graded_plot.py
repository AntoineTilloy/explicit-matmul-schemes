#!/usr/bin/env python3
"""Three-range public rational-scheme figure.

Render on Atlas. Reads the public catalogue; writes exponents.svg/png/pdf
and a small provenance manifest. Does not read or alter characteristic-2 data.
No data-point interpolation. Arbitrarily large n and ranks stay Python ints;
only their logarithms are converted to plotting coordinates.

Integration: draw_catalogue(catalogue_path, figures_directory).
CLI: python graded_plot.py --catalogue catalogue.json --out figures
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

LITA32_RANK = 14196
LITA32_EXPONENT = math.log(LITA32_RANK)/math.log(32)
MEDIUM_LIMIT = 30000000000
COLORS = dict(cw='#176c99', z4='#7b4ea3', small='#d9822b',
              rectangular='#348568', graded='#c53e55', reference='#7c8794')
STYLES = {
    'cw': ('CW constructions', 'o', 27),
    'z4': ('CW quotient constructions', 's', 28),
    'small': ('Small rational schemes', 'D', 32),
    'rectangular': ('Rectangular compositions / LITA refresh', '^', 35),
    'graded': ('Graded constructions', '*', 125),
}


def log10_integer(value):
    """Avoid float(value), which overflows for the largest finite schemes."""
    if value <= 0:
        raise ValueError('Logarithms require positive integers')
    shift = max(0, value.bit_length()-53)
    return math.log10(value >> shift) + shift*math.log10(2)


def approximate_size(value):
    if value < 1000000000000:
        return f'{value:,}'
    exponent = int(math.floor(log10_integer(value)))
    mantissa = 10**(log10_integer(value)-exponent)
    if mantissa >= 9.995:
        mantissa /= 10
        exponent += 1
    return rf'${mantissa:.2f}\times10^{{{exponent}}}$'


def classify(entry, n):
    provenance = json.dumps(entry.get('provenance', {}), sort_keys=True).lower()
    if 'graded' in provenance:
        return 'graded'
    if n <= 32:
        return 'small'
    if 'z4-quotient' in provenance:
        return 'z4'
    if 'rectangular' in provenance or 'upstream primitive refresh' in provenance:
        return 'rectangular'
    return 'cw'


def prepare_points(entries):
    points = []
    for entry in entries:
        if entry.get('field') != 'Q':
            continue
        n, rank = int(entry['n']), int(entry['rank'])
        if n <= 1 or rank <= 0:
            raise ValueError('Catalogue plot needs n>1 and positive rank')
        x = log10_integer(n)
        family = classify(entry, n)
        if family == 'graded' and not (entry.get('root_object') and entry.get('verification')):
            raise ValueError('Graded plot points require a catalogue root and verification description')
        points.append(dict(n=n, rank=rank, x=x, exponent=log10_integer(rank)/x,
                           family=family, root=entry.get('root_object'),
                           sub265=(rank**20 < n**53) if family == 'graded' else False))
    points.sort(key=lambda p: p['n'])
    if len({p['n'] for p in points}) != len(points):
        raise ValueError('Expected one best rational catalogue entry per size')
    if not any(p['n'] == 32 and p['rank'] == LITA32_RANK for p in points):
        raise ValueError('Missing pinned n=32, rank=14196 rational reference')
    return points


def exponent_label(point):
    # Cropping can make the strict inequality invisible at six decimal places.
    if point['sub265'] and point['exponent'] >= 2.6499995:
        return r'$\tau<2.65$ (exact integer check)'
    return f"exponent {point['exponent']:.6f}"


def annotate(ax, point, text, location, color='#354553'):
    ax.annotate(text, (point['x'], point['exponent']), xytext=location,
                textcoords='axes fraction', fontsize=8.8, color=color,
                ha='left', va='bottom',
                bbox=dict(boxstyle='round,pad=.23', fc='white', ec='none', alpha=.92),
                arrowprops=dict(arrowstyle='-', color='#9da8b0', lw=.8),
                annotation_clip=True)


def wide_ticks(upper):
    raw = max(1.0, upper/6)
    magnitude = 10**math.floor(math.log10(raw))
    step = next(v*magnitude for v in (1, 2, 5, 10) if v*magnitude >= raw)
    ticks = [1]
    ticks.extend(range(int(step), int(upper)+1, int(step)))
    return sorted(set(ticks))


def draw_catalogue(catalogue_path, output_directory):
    if hasattr(sys, 'set_int_max_str_digits'):
        sys.set_int_max_str_digits(0)
    catalogue_path, output_directory = Path(catalogue_path), Path(output_directory)
    raw = catalogue_path.read_bytes()
    points = prepare_points(json.loads(raw))
    old = [p for p in points if p['family'] != 'graded']
    graded = [p for p in points if p['family'] == 'graded']
    small = [p for p in points if p['n'] <= 16384]
    medium = [p for p in points if p['n'] <= MEDIUM_LIMIT]
    if not small or not medium or not old:
        raise ValueError('The small and medium panels require rational entries')

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FuncFormatter
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'svg.fonttype': 'none', 'axes.spines.top': False,
                         'axes.spines.right': False, 'pdf.fonttype': 42})
    fig, axes = plt.subplots(1, 3, figsize=(18, 7.4),
                             gridspec_kw={'width_ratios': [1.08, 1.0, 1.08]})
    full_max = max(p['x'] for p in points)
    upper = max(math.log10(MEDIUM_LIMIT), full_max+max(.5, .045*full_max))
    views = [(small, 'Small sizes', (.95, math.log10(20000)),
              (2.715, 2.802), [1, 2, 3, 4]),
             (medium, 'Medium sizes', (.95, math.log10(MEDIUM_LIMIT)),
              (2.685, 2.802), [1, 3, 5, 7, 9]),
             (points, 'Large sizes', (.95, upper),
              (min(2.635, min(p['exponent'] for p in points)-.012), 2.802), wide_ticks(upper))]
    for ax, (subset, title, xlim, ylim, ticks) in zip(axes, views):
        for family, (label, marker, size) in STYLES.items():
            group = [p for p in subset if p['family'] == family]
            if group:
                ax.scatter([p['x'] for p in group], [p['exponent'] for p in group],
                           c=COLORS[family], marker=marker, s=size,
                           linewidths=.5 if family == 'graded' else 0,
                           edgecolors='white' if family == 'graded' else 'none',
                           zorder=6 if family == 'graded' else 4)
        ax.axhline(LITA32_EXPONENT, color=COLORS['reference'], ls='--', lw=1, zorder=2)
        ax.set_xlim(*xlim)
        # Preserve the old detail while allowing any future point below it.
        lower = min(ylim[0], min(p['exponent'] for p in subset)-.008)
        ax.set_ylim(lower, ylim[1])
        ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: rf'$10^{{{int(x)}}}$'))
        ax.grid(axis='y', color='#e9edf1', lw=.7)
        ax.set_axisbelow(True)
        ax.set_title(title, loc='left', fontsize=12, pad=13)
        ax.set_xlabel(r'Square matrix size $n$ (logarithmic scale)', labelpad=10)
    axes[0].set_ylabel(r'Effective exponent $\log_n R(n)$ — lower is better', labelpad=10)
    axes[1].axhline(2.7, color='#479b7d', ls=':', lw=.9, zorder=2)
    axes[1].text(.02, (2.7-axes[1].get_ylim()[0])/(axes[1].get_ylim()[1]-axes[1].get_ylim()[0])+.015,
                 '2.7', transform=axes[1].transAxes, color='#38856c', fontsize=8.5)
    axes[2].axhline(2.65, color='#ba4053', ls=':', lw=1.1, zorder=2)
    axes[2].text(.02, (2.65-axes[2].get_ylim()[0])/(axes[2].get_ylim()[1]-axes[2].get_ylim()[0])+.015,
                 '2.65 target', transform=axes[2].transAxes, color='#ba4053', fontsize=8.8)

    for point in small:
        if point['n'] <= 32:
            axes[0].annotate(f"{point['n']}: {point['rank']:,}",
                             (point['x'], point['exponent']), xytext=(7, 4),
                             textcoords='offset points', color='#8a4f13', fontsize=8.3)
    n1024 = next((p for p in small if p['n'] == 1024), None)
    if n1024:
        annotate(axes[0], n1024,
                 f"1,024: {n1024['rank']:,}\n{exponent_label(n1024)}", (.47, .76), COLORS['z4'])
    best_small = min(small, key=lambda p: p['exponent'])
    annotate(axes[0], best_small,
             f"{best_small['exponent']:.6f} at n = {best_small['n']:,}", (.40, .09))
    best_old = min((p for p in medium if p['family'] != 'graded'), key=lambda p: p['exponent'])
    annotate(axes[1], best_old,
             f"Best exponent: {best_old['exponent']:.6f}\nn = {best_old['n']:,}", (.19, .20))
    if graded:
        best_graded = min(graded, key=lambda p: p['exponent'])
        annotate(axes[2], best_graded,
                 f"Graded construction\n{exponent_label(best_graded)}\nn ≈ {approximate_size(best_graded['n'])}",
                 (.43, .30), COLORS['graded'])

    legend = [Line2D([], [], ls='none', marker=marker, markersize=8 if family == 'graded' else 6,
                     color=COLORS[family], label=label)
              for family, (label, marker, _) in STYLES.items()
              if any(p['family'] == family for p in points)]
    legend.append(Line2D([], [], ls='--', lw=1, color=COLORS['reference'],
                         label=f'32-LITA reference: {LITA32_EXPONENT:.6f}'))
    fig.legend(handles=legend, loc='lower center', bbox_to_anchor=(.52, .135),
               ncol=3, frameon=False, fontsize=9.3, columnspacing=2.0, handletextpad=.65)
    fig.suptitle('Explicit square matrix multiplication schemes', x=.055, y=.965,
                 ha='left', fontsize=20, weight='bold')
    fig.text(.055, .905,
             f'{len(points)} saved rational schemes over R and C  •  {len(graded)} graded constructions  •  One best saved scheme per size',
             fontsize=11, color='#526171')
    fig.text(.055, .072, 'Each marker is a saved finite construction. Lines between construction points are not drawn.',
             fontsize=9.5, color='#526171')
    fig.text(.055, .039, '32-LITA attains its reference exponent at powers of 32; other finite sizes incur blocking and padding costs.',
             fontsize=9.5, color='#526171')
    fig.subplots_adjust(left=.055, right=.985, top=.805, bottom=.29, wspace=.25)
    output_directory.mkdir(parents=True, exist_ok=True)
    for extension in ('svg', 'png', 'pdf'):
        fig.savefig(output_directory/f'exponents.{extension}', dpi=180, facecolor='white')
    plt.close(fig)
    manifest = dict(catalogue_sha256=hashlib.sha256(raw).hexdigest(),
                    source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    plotted_field='Q (valid over R and C)', points=len(points),
                    graded_points=len(graded), x_coordinate='log10(n), computed from integer leading bits',
                    data_interpolation=False, lita32_rank=str(LITA32_RANK),
                    lita32_exponent=LITA32_EXPONENT,
                    best_earlier=dict(n=str(best_old['n']), rank=str(best_old['rank']), root=best_old['root']),
                    graded=[dict(n=str(p['n']), rank=str(p['rank']), root=p['root'],
                                 exponent_for_display=p['exponent'], strict_sub265=p['sub265']) for p in graded])
    (output_directory/'exponents.metadata.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--catalogue', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if hasattr(sys, 'set_int_max_str_digits'):
        sys.set_int_max_str_digits(0)
    print(json.dumps(draw_catalogue(args.catalogue, args.out), indent=2))


if __name__ == '__main__':
    main()
