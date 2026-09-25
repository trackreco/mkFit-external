#!/usr/bin/env python3
"""Turn `seedfind --truth` tables into ROOT graphs versus |eta|.

    truth2root.py OUT.root  0.2=pt0.2.txt  0.5=pt0.5.txt  ...

For each table (label = its pT threshold in GeV) it writes, as TGraphAsymmErrors
with Clopper-Pearson 68 % intervals where the quantity is a fraction:

    eff_<lbl>     found / findable                     vs the sim track's |eta|
    fake_<lbl>    fake / (true + fake), decidable only  vs the quad's |eta|
    undec_<lbl>   undecidable / all quads               vs the quad's |eta|
    dup_<lbl>     extra true quads per found track      vs the sim track's |eta|
                  (error sqrt(extra) / found)
    dupf_<lbl>    found tracks with >= 2 true quads / found

and one TMultiGraph per quantity (mg_eff, mg_fake, mg_undec, mg_dup, mg_dupf)
holding every threshold, coloured 1 / 4 / 2 / 8 in the order given, so a JSROOT
page can draw one object per panel. Bins with an empty denominator are left out.
"""
import sys
import ROOT

ROOT.gROOT.SetBatch(True)
COLORS = [1, 4, 2, 8, 6, 9]
MARKERS = [20, 21, 22, 23, 33, 34]


def read(fn):
    rows = []
    for line in open(fn):
        if line.startswith('#') or not line.strip():
            continue
        rows.append([float(x) for x in line.split()])
    return rows


def frac_graph(name, rows, num, den):
    g = ROOT.TGraphAsymmErrors()
    g.SetName(name)
    for r in rows:
        lo, hi, n, d = r[0], r[1], r[num], r[den] if isinstance(den, int) else den(r)
        if d <= 0:
            continue
        e = n / d
        elo = e - ROOT.TEfficiency.ClopperPearson(int(d), int(n), 0.683, False)
        ehi = ROOT.TEfficiency.ClopperPearson(int(d), int(n), 0.683, True) - e
        i = g.GetN()
        g.SetPoint(i, 0.5 * (lo + hi), e)
        g.SetPointError(i, 0.5 * (hi - lo), 0.5 * (hi - lo), elo, ehi)
    return g


def dup_graph(name, rows):
    g = ROOT.TGraphAsymmErrors()
    g.SetName(name)
    for r in rows:
        lo, hi, found, extra = r[0], r[1], r[3], r[4]
        if found <= 0:
            continue
        v, e = extra / found, (extra ** 0.5) / found
        i = g.GetN()
        g.SetPoint(i, 0.5 * (lo + hi), v)
        g.SetPointError(i, 0.5 * (hi - lo), 0.5 * (hi - lo), min(e, v), e)
    return g


def main():
    out = ROOT.TFile(sys.argv[1], 'RECREATE')
    titles = {
        'eff': ';|#eta| of the sim track;seeding efficiency',
        'fake': ';|#eta| of the quadruplet;fake rate (decidable quads)',
        'undec': ';|#eta| of the quadruplet;undecidable fraction',
        'dup': ';|#eta| of the sim track;extra true quads per found track',
        'dupf': ';|#eta| of the sim track;found tracks with #geq 2 true quads',
    }
    mgs = {k: ROOT.TMultiGraph('mg_' + k, t) for k, t in titles.items()}
    for i, arg in enumerate(sys.argv[2:]):
        lbl, fn = arg.split('=', 1)
        rows = read(fn)
        tag = lbl.replace('.', 'p')
        gs = {
            'eff': frac_graph('eff_' + tag, rows, 3, 2),
            'fake': frac_graph('fake_' + tag, rows, 8, lambda r: r[7] + r[8]),
            'undec': frac_graph('undec_' + tag, rows, 9, 6),
            'dup': dup_graph('dup_' + tag, rows),
            'dupf': frac_graph('dupf_' + tag, rows, 5, 3),
        }
        for k, g in gs.items():
            g.SetTitle('pT > %s GeV' % lbl)
            g.SetLineColor(COLORS[i])
            g.SetMarkerColor(COLORS[i])
            g.SetMarkerStyle(MARKERS[i])
            g.SetMarkerSize(1.0)
            g.Write()
            mgs[k].Add(g.Clone(), 'P')
    for mg in mgs.values():
        mg.Write()
    out.Close()


if __name__ == '__main__':
    main()
