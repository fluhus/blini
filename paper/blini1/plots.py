import csv
from collections import defaultdict

import numpy as np
from matplotlib import pyplot as plt
from myplot import ctx

OUT_DIR = '../testdata/output'
SIZE_RATIO = [0.5, 0.75]


def load_table(f: str):
    d = defaultdict(list)
    for row in csv.DictReader(open(f)):
        for k, v in row.items():
            d[k].append(v)
    for k in d:
        try:
            d[k] = [float(x) for x in d[k]]
        except ValueError:
            continue
    if 'time_seconds_1' in d:
        d['time_seconds'] = list(zip(*[d[f'time_seconds_{i}'] for i in range(1, 6)]))
        d['time_seconds'] = [np.mean(x) for x in d['time_seconds']]
    return d


def plot_search_small():
    d = load_table('results/results_search.txt')
    d['name'] = [x.replace(' ', '\n') for x in d['name']]
    with ctx(f'{OUT_DIR}/search_time', sizeratio=SIZE_RATIO):
        plt.bar(d['name'][:2], d['time_seconds'][:2])
        plt.bar(d['name'][2:4], d['time_seconds'][2:4])
        plt.bar(d['name'][4:], d['time_seconds'][4:])
        plt.yscale('log')
        plt.ylabel('Average time (s)')
        plt.xticks(rotation=90)
    with ctx(f'{OUT_DIR}/search_found', sizeratio=SIZE_RATIO):
        plt.bar(d['name'][:2], d['source_found'][:2])
        plt.bar(d['name'][2:4], d['source_found'][2:4])
        plt.bar(d['name'][4:], d['source_found'][4:])
        plt.ylabel('Successful matches (out of 100)')
        plt.xticks(rotation=90)
    with ctx(f'{OUT_DIR}/search_others', sizeratio=SIZE_RATIO):
        plt.bar(d['name'][:2], d['others_found'][:2])
        plt.bar(d['name'][2:4], d['others_found'][2:4])
        plt.bar(d['name'][4:], d['others_found'][4:])
        plt.ylabel('Non-source matches ')
        plt.xticks(rotation=90)


def plot_search_big():
    d = load_table('results/results_search_big.txt')
    d['name'] = [x.replace(' (', '\n(') for x in d['name']]
    with ctx('results/search_big', sizeratio=0.75):
        plt.bar(d['name'][:3], d['time_seconds'][:3])
        plt.bar(d['name'][3:5], d['time_seconds'][3:5])
        plt.bar(d['name'][5:], d['time_seconds'][5:])
        plt.text(5, d['time_seconds'][5] * 1.02, 'X', ha='center')
        plt.yscale('log')
        plt.ylabel('Time (s)')
        plt.xticks(rotation=90)


def plot_clust(tag: str):
    d = load_table(f'results/results_clust_{tag}.txt')
    d['name'] = [x.replace(' (', '\n(') for x in d['name']]
    with ctx(f'{OUT_DIR}/clust_{tag}_time', sizeratio=SIZE_RATIO):
        plt.bar(d['name'][:4], d['time_seconds'][:4])
        plt.bar(d['name'][4:], d['time_seconds'][4:])
        plt.ylabel('Average time (s)')
        plt.xticks(rotation=90)
    with ctx(f'{OUT_DIR}/clust_{tag}_mem', sizeratio=SIZE_RATIO):
        plt.bar(d['name'][:4], d['max_mem_mb'][:4])
        plt.bar(d['name'][4:], d['max_mem_mb'][4:])
        plt.ylabel('Max memory (MB)')
        plt.xticks(rotation=90)
    with ctx(f'{OUT_DIR}/clust_{tag}_ari', sizeratio=SIZE_RATIO):
        plt.bar(d['name'][:4], d['ari'][:4])
        plt.bar(d['name'][4:], d['ari'][4:])
        plt.ylabel('Adjusted Rand-Index')
        ylim = [min(d['ari']), max(d['ari'])]
        stretch = (ylim[1] - ylim[0]) * 0.2
        ylim[0] -= stretch
        ylim[1] += stretch
        plt.ylim(ylim)
        plt.xticks(rotation=90)
    with ctx(f'{OUT_DIR}/clust_{tag}_nclust', sizeratio=SIZE_RATIO):
        plt.bar(d['name'][:4], d['n_clusters'][:4])
        plt.bar(d['name'][4:], d['n_clusters'][4:])
        plt.ylabel('Number of clusters')
        plt.xticks(rotation=90)


plt.style.use('ggplot')
plot_search_small()
plot_search_big()
plot_clust('frag')
plot_clust('snps')
