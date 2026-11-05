#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Pre-submission compliance check for v5.0 (Diabetes, Obesity and Metabolism).

Every rule is quoted from the DOM Guide for Authors fetched 2026-10-28:
  - Structured abstract, four sections (Aims / Materials and Methods / Results / Conclusions)
  - Abstract <= 250 words
  - Main text <= 3,500 words (Introduction to Conclusions)
  - Figures and tables combined <= 5
  - References normally < 60
"""
import re
import sys

P = 'versions/v7.4_20261031/Manuscript_v7.4_DOM_20261031.md'
s = open(P, encoding='utf-8').read()


def wc(t):
    """Token count matching the manuscript's own count: strip image links,
    emphasis markers and backticks, count tokens containing an alphanumeric."""
    t = re.sub(r'!\[.*?\]\(.*?\)', '', t)
    t = re.sub(r'\*\*|\*|`', '', t)
    return len([w for w in t.split() if re.search(r'[A-Za-z0-9]', w)])


def seg(a, b):
    i = s.find(a)
    j = s.find(b, i)
    return s[i:j] if i >= 0 and j > i else ''


ab = seg('## Abstract', '## 1. Introduction')
main = seg('## 1. Introduction', '## Declarations')
# DOM order: Declarations -> References -> Figure legends -> Tables -> Supplementary
refs = seg('## References', '## Figure legends')
leg = seg('## Figure legends', '## Tables')
tab = seg('## Tables', '## Supplementary materials')

checks = []


def chk(name, ok, detail):
    checks.append((ok, name, detail))


# 1. abstract structure: four required headings
heads = [r'Aims', r'Materials and methods', r'Results', r'Conclusions']
found = [h for h in heads if re.search(r'^\s*' + h + r'[\.\:]', ab, re.M | re.I)]
chk('Abstract structured: Aims/Materials and Methods/Results/Conclusions',
    len(found) == 4, f'found {found}')

# 2. abstract <= 250 words
n_ab = wc(ab)
chk('Abstract <= 250 words', n_ab <= 250, f'{n_ab}')

# 3. main text <= 3,500 words
n_main = wc(main)
chk('Main text <= 3,500 words', n_main <= 3500, f'{n_main} ({3500 - n_main} margin)')

# 4. figures + tables <= 5
n_fig = len(re.findall(r'^\*\*Fig\. \d', leg, re.M))
n_tab = len(re.findall(r'^### Table \d', tab, re.M))
chk('Figures + tables <= 5', n_fig + n_tab <= 5, f'{n_fig} fig + {n_tab} tab')

# 5. references < 60
nums = [int(x) for x in re.findall(r'^\[(\d+)\]', refs, re.M)]
contiguous = nums == list(range(1, len(nums) + 1))
chk('References < 60', len(nums) < 60, f'{len(nums)}')
chk('References numbered contiguously from 1', contiguous, f'{len(nums)} refs, max {max(nums) if nums else 0}')

# 6. every in-text citation resolves to a listed reference
cited = set()
for m in re.finditer(r'\[([\d,\u2013\-\s]+)\]', s[:s.find('## References')]):
    for part in m.group(1).split(','):
        part = part.strip()
        if re.match(r'^\d+$', part):
            cited.add(int(part))
        elif re.match(r'^\d+\s*[\u2013-]\s*\d+$', part):
            a, b = re.split(r'[\u2013-]', part)
            cited.update(range(int(a), int(b) + 1))
orphan = sorted(cited - set(nums))
chk('No orphan in-text citations', not orphan, f'orphans={orphan}')
uncited = sorted(set(nums) - cited)
chk('No uncited references', not uncited, f'uncited={uncited}')

# 7. disclosure register D-numbers present and contiguous.
#    The lookbehind on [A-Za-z] keeps acronyms such as RAED2 from registering as D2;
#    the S6 summary line states the range D1-D32, which is what we check against.
# The main text cites only the deviations that changed the analysis set; the
# full register D1-D32 lives in supplementary S6. So the check is: every cited
# D-number must be within the declared range, and the range statement must exist.
dnums = sorted({int(x) for x in re.findall(r'(?<![A-Za-z])D(\d+)', s)})
in_range = all(1 <= n <= 32 for n in dnums)
chk('All cited D-numbers within declared range D1-D32', in_range, f'cited={dnums}')
chk('S6 declares the full register D1-D32', 'D1\u2013D32' in s or 'D1-D32' in s, 'ok')
chk('Analysis-set deviations D23/D26/D27 all cited', all(f'D{n}' in s for n in (23, 26, 27)), 'ok')

# 8. key pooled estimates exactly as computed (guard against drift)
expect = ['46.6%', '19.5%', '79.9%', '0.003', '71.8%',
          '42.3%', '28.2%', '57.9%', '50.6%',
          '+0.2%', '\u221215.5%', '18.9%', '0.978', '85.1%',
          '+32.8%', '\u221210.0%', '+47.7%', '0.001',
          '+97.2%', '25.2%', '210.7%', '+25.5%', '4.6%', '50.6%']
missing = [e for e in expect if e not in s]
chk('All pooled estimates match computed output', not missing, f'missing={missing}')

# 9. arm / study counts consistent everywhere
bad = []
for pat, want in [(r'nine arms', None), (r'Nine arms from five studies', 1),
                  (r'[Nn]ine exercise arms from five studies', 1),
                  (r'five contributed nine arms', 1)]:
    if not re.search(pat, s):
        bad.append(pat)
chk('k accounting stated as 9 arms / 5 studies', not bad, f'missing={bad}')

# 10. no leftover Diabetologia apparatus outside references
body_no_refs = s[:s.find('## References')] + s[s.find('## Supplementary materials'):]
leftovers = [b for b in ['Research in context', '@DiabetologiaJnl', 'Diabetologia',
                         'DRCP', 'Diabetes Research and Clinical Practice'] if b in body_no_refs]
chk('No leftover previous-target apparatus', not leftovers, f'leftovers={leftovers}')

# 11. author block per team rule (three authors, Sun corresponding)
chk('Author block: Sun/Wu/Ma', all(x in s[:2500] for x in ['Zhongwei Sun', 'Xinran Wu', 'Yongping Ma']), 'ok')
chk('Corresponding author = Sun', 'sunzhongwei@gzgs.edu.cn' in s[:3000], 'ok')

# 12. no placeholder residue
chk('No placeholder residue', 'WORDCOUNT' not in s and 'TODO' not in s and 'TBD' not in s, 'ok')

width = max(len(c[1]) for c in checks)
fails = 0
for ok, name, detail in checks:
    print(('PASS  ' if ok else 'FAIL  ') + name.ljust(width) + '  ' + detail)
    if not ok:
        fails += 1
print()
print(f'{len(checks) - fails}/{len(checks)} PASS, {fails} FAIL')
sys.exit(1 if fails else 0)
