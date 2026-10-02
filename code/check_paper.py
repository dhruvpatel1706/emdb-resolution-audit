"""Pre-submission checks for paper/main.tex. Exit status 1 if any check fails.

1. numbers.tex, the tables and the figures are regenerated from the saved results, and numbers.tex must not change.
2. Every macro the text uses for a number is defined in numbers.tex.
3. No literal decimal number in the prose except a fixed list of definitions (thresholds, bands, software versions).
4. Style: no em-dashes, no first-person plural, none of the phrases a referee flagged in the rejected version, none of
   the usual filler words.
5. No placeholder text (blocks staging until the repository link is filled in).
6. The PDF is newer than every file it is built from.

    python3 code/check_paper.py
"""
import hashlib, os, re, subprocess, sys

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
P = os.path.join(ROOT, "paper")
fails = []


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


GEN = ["numbers.tex", "tab_selection.tex", "tab_coverage.tex", "tab_agreement.tex", "tab_band_method.tex", "tab_unindexed.tex",
       "fig1_ratio_distributions.pdf", "fig2_stated_vs_corrected.pdf", "fig3_tail_by_band_and_method.pdf"]
before = {f: sha(os.path.join(P, f)) for f in GEN}
subprocess.run([sys.executable, os.path.join(ROOT, "code", "make_paper_inputs.py"),
                os.path.join(ROOT, "data", "raw", "emdb_index_2026-10-01.csv")], check=True, capture_output=True)
changed = [f for f in GEN if sha(os.path.join(P, f)) != before[f]]
if changed:
    fails.append(f"regenerated files differ from the ones in use: {changed}")
for cmd in (["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "main.tex"], ["bibtex", "main"],
            ["pdflatex", "-interaction=nonstopmode", "main.tex"], ["pdflatex", "-interaction=nonstopmode", "main.tex"]):
    if subprocess.run(cmd, cwd=P, capture_output=True).returncode != 0:
        fails.append(f"build step failed: {' '.join(cmd)}")
        break
log = open(os.path.join(P, "main.log"), encoding="latin-1").read()
for pat in ("Overfull \\hbox", "undefined", "\n! "):
    if pat in log:
        fails.append(f"LaTeX log contains {pat.strip()!r}")

tex = open(os.path.join(P, "main.tex"), encoding="utf-8").read()
body = tex[tex.index("\\begin{document}"):tex.index("\\bibliographystyle")]
defined = set(re.findall(r"\\newcommand\{\\(\w+)\}", open(os.path.join(P, "numbers.tex")).read()))
known_latex = {"maketitle", "begin", "end", "noindent", "textbf", "section", "subsection", "citep", "citet", "label", "ref",
               "url", "texttt", "small", "footnotesize", "ttfamily", "normalfont", "centering", "caption", "toprule",
               "midrule", "bottomrule", "cmidrule", "multicolumn", "tabinput", "includegraphics", "linewidth", "AA",
               "setlength", "tabcolsep", "raggedright", "arraybackslash", "item", "itemsep", "emph", "bibliographystyle"}
used = set(re.findall(r"\\([A-Za-z]+)", body)) - known_latex
missing = sorted(u for u in used if u not in defined and u[0].islower() and any(c.isupper() for c in u[1:]))
if missing:
    fails.append(f"macros used but not defined: {missing}")

prose = re.sub(r"\\(texttt|url|label|ref|includegraphics|tabinput|cite[pt])(\[[^]]*\])?\{[^}]*\}", " ", body)
prose = re.sub(r"\\begin\{tabular\}.*?\\end\{tabular\}", " ", prose, flags=re.S)
allowed = {"0.143", "1.10", "1.25", "1.50", "2.00", "1.05", "0.01", "3.14", "2.4", "3.11"}
lits = sorted({m for m in re.findall(r"(?<![\w.\\-])\d+\.\d+", prose) if m not in allowed})
if lits:
    fails.append(f"literal decimals in prose (should be macros): {lits}")

if "\u2014" in tex or "---" in tex:
    fails.append("em-dash present")
words = re.findall(r"[A-Za-z']+", prose.lower())
for w in ("we", "our", "us", "ours"):
    if w in words:
        fails.append(f"first-person plural '{w}' in a single-author paper")
for phrase in ("spine", "absorb", " lags", "bounded delta", "ratif", "weak positive control", "runs on a cpu",
               "delve", "leverage", "utilize", "robust", "crucial", "pivotal", "notably", "moreover", "furthermore",
               "comprehensive", "landscape", "underscore", "it is worth noting", "plays a key role", "sheds light"):
    if phrase in prose.lower():
        fails.append(f"phrase '{phrase.strip()}' in text")

if re.search(r"TO BE ADDED|TODO|XXX|PLACEHOLDER", tex):
    fails.append("placeholder text present (repository link not filled in)")

pdf = os.path.join(P, "main.pdf")
srcs = [os.path.join(P, f) for f in os.listdir(P) if f.endswith((".tex", ".bib", ".pdf")) and f != "main.pdf"]
if not os.path.exists(pdf) or any(os.path.getmtime(s) > os.path.getmtime(pdf) for s in srcs):
    fails.append("main.pdf is missing or older than its sources (rebuild)")

print(f"macros used: {len(used & defined)} of {len(defined)} defined")
for f in fails:
    print("FAIL:", f)
print("RESULT:", "FAIL" if fails else "PASS")
sys.exit(1 if fails else 0)
