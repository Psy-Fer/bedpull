#!/usr/bin/env python3
"""Score bedpull vs liftOver extraction accuracy against the GIAB SV truth set.

Each tool may produce more than one results.tsv (one per HG002 haplotype —
see run_bedpull_pipeline.sh / run_liftover_pipeline.sh) and, for bedpull,
more than one hit within a single haplotype run. For each truth window we
take whichever hit (across all haplotype files) comes closest to the
expected extracted length, since a het SV is only real on one haplotype and
we don't try to map VCF GT phase to assembly haplotype identity upstream.

Expected extracted length is derived directly from the window size and the
REF/ALT length delta recorded in truth.tsv:

    expected_length = window_len + (alt_len - ref_len)

This is positive for insertions (assembly has extra bases the reference
window doesn't) and negative-adjusting for deletions (assembly is missing
bases the reference window spans) — one formula, no SVTYPE branching.

Outputs a markdown summary report (overall + by SVTYPE + by size bucket +
worst-offenders) and a full per-window TSV for further analysis.

Precision/F1/ROC (optional, needs --neg-truth + --neg-*-results): without a
negative-control set, "precision" is meaningless here — every window in
truth.tsv is a window we already know has a real SV, so a tool can never be
observed producing a spurious result on it. select_negative_control_windows.py
generates windows matched in size but far from any known SV; scoring those
the same way gives real TP/FP/FN/TN counts:

    TP = true-SV window,     extracted within tolerance of the true length
    FN = true-SV window,     missing or outside tolerance (unmapped counts here)
    FP = negative window,    extracted something outside tolerance of window_len
    TN = negative window,    extracted within tolerance of window_len

Negative windows with NO output at all are excluded from the confusion
matrix entirely (reported separately) — failing to extract anything isn't
the same claim as extracting the wrong length, so it shouldn't count as
either a false or a true negative.
"""

import argparse
import csv
import subprocess
from collections import defaultdict

SIZE_BUCKETS = [
    (50, 99, "50-99bp"),
    (100, 499, "100-499bp"),
    (500, 4999, "500-4999bp"),
    (5000, float("inf"), ">=5000bp"),
]


def size_bucket(svlen):
    for lo, hi, label in SIZE_BUCKETS:
        if lo <= svlen <= hi:
            return label
    return "other"


def load_truth(path):
    truth = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            row["svlen"] = int(row["svlen"])
            row["ref_len"] = int(row["ref_len"])
            row["alt_len"] = int(row["alt_len"])
            row["window_start"] = int(row["window_start"])
            row["window_end"] = int(row["window_end"])
            row["window_len"] = row["window_end"] - row["window_start"]
            row["expected_length"] = row["window_len"] + (row["alt_len"] - row["ref_len"])
            truth[row["name"]] = row
    return truth


def load_results(paths):
    """Return {name: [row, ...]} merged across all given results.tsv files."""
    by_name = defaultdict(list)
    for path in paths:
        with open(path) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                row["extracted_len"] = int(row["extracted_len"]) if row["extracted_len"] else 0
                by_name[row["name"]].append(row)
    return by_name


def best_hit(rows, expected_length):
    """Pick the 'ok' row closest to expected_length; None if all unmapped."""
    ok_rows = [r for r in rows if r["status"] == "ok"]
    if not ok_rows:
        return None
    return min(ok_rows, key=lambda r: abs(r["extracted_len"] - expected_length))


def score_tool(truth, results_by_name, tolerance_frac, tolerance_floor):
    """Return a list of per-window score dicts for one tool."""
    scored = []
    for name, t in truth.items():
        rows = results_by_name.get(name, [])
        expected = t["expected_length"]
        tolerance = max(tolerance_floor, tolerance_frac * expected)
        hit = best_hit(rows, expected)
        if hit is None:
            unmapped_reasons = sorted({r["unmapped_reason"] for r in rows if r["unmapped_reason"]})
            scored.append(
                {
                    "name": name,
                    "svtype": t["svtype"],
                    "svlen": t["svlen"],
                    "size_bucket": size_bucket(t["svlen"]),
                    "expected_length": expected,
                    "extracted_length": 0,
                    "abs_error": expected,
                    "status": "unmapped",
                    "pass": False,
                    "looks_like_no_variant": False,
                    "unmapped_reason": "; ".join(unmapped_reasons) if unmapped_reasons else "no hit",
                }
            )
            continue
        abs_error = abs(hit["extracted_len"] - expected)
        passed = abs_error <= tolerance
        # A failed extraction that instead matches the *unvaried* window length
        # looks like "no variant found here" rather than "wrong-sized variant" —
        # worth distinguishing when the truth set and the extraction target
        # might be different genomes (see the confirmation-leg discussion in
        # README.md): that pattern is often the truth set expecting a variant
        # that simply isn't present in the sequence actually being extracted,
        # not a tool failure.
        looks_like_no_variant = False
        if not passed:
            no_variant_tolerance = max(tolerance_floor, tolerance_frac * t["window_len"])
            looks_like_no_variant = abs(hit["extracted_len"] - t["window_len"]) <= no_variant_tolerance
        scored.append(
            {
                "name": name,
                "svtype": t["svtype"],
                "svlen": t["svlen"],
                "size_bucket": size_bucket(t["svlen"]),
                "expected_length": expected,
                "extracted_length": hit["extracted_len"],
                "abs_error": abs_error,
                "status": "ok",
                "pass": passed,
                "looks_like_no_variant": looks_like_no_variant,
                "unmapped_reason": "",
            }
        )
    return scored


def confusion_matrix(pos_scored, neg_scored):
    """TP/FN from the positive set, FP/TN from the negative set, at the
    tolerance already baked into each entry's 'pass' field. Unmapped
    negative-control windows are excluded (see module docstring)."""
    tp = sum(1 for r in pos_scored if r["pass"])
    fn = len(pos_scored) - tp
    neg_considered = [r for r in neg_scored if r["status"] == "ok"]
    fp = sum(1 for r in neg_considered if not r["pass"])
    tn = len(neg_considered) - fp
    n_neg_excluded = len(neg_scored) - len(neg_considered)

    precision = tp / (tp + fp) if (tp + fp) else float("nan")
    recall = tp / (tp + fn) if (tp + fn) else float("nan")
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else float("nan")
    return {
        "tp": tp,
        "fn": fn,
        "fp": fp,
        "tn": tn,
        "n_neg_excluded_unmapped": n_neg_excluded,
        "precision": precision,
        "recall": recall,
        "f1": f1,
    }


def detection_score(row, truth_row):
    """Deviation of extracted length from the plain reference-window length
    (window_end - window_start) — deliberately *not* the class-specific
    expected_length. expected_length differs in meaning between a true-SV
    window (matching it means "got the SV right") and a negative-control
    window (matching it means "correctly found nothing"), so sweeping a
    threshold on relative-to-expected error doesn't trace a real ROC: at
    any threshold, "small error" means opposite things for the two classes,
    which is why an earlier version of this function produced a nonsense
    curve (FPR was already >85% at the strictest threshold).

    Deviation from window_len, by contrast, is one detection question asked
    the same way of both classes: "does the extracted length look different
    from a plain, no-indel reference span?" A true SV should trip this; a
    negative control shouldn't. This says nothing about whether a detected
    SV's *size* is correct — that's what confusion_matrix()/F1 are for.

    None for unmapped rows — "no output" never counts as a detection, at
    any threshold.
    """
    if row["status"] != "ok":
        return None
    window_len = truth_row["window_end"] - truth_row["window_start"]
    return abs(row["extracted_length"] - window_len) / max(window_len, 1)


def roc_points(pos_scored, pos_truth, neg_scored, neg_truth):
    """Empirical detection ROC: TPR/FPR, as a function of how large a
    deviation from the plain reference-window length counts as "detected",
    at every distinct deviation value observed across both sets (the
    standard discrete-scores ROC construction, as in
    sklearn.metrics.roc_curve). Unmapped negative windows are excluded from
    the FPR denominator for the same reason they're excluded from the
    confusion matrix — see the module docstring.

    A window is "detected" (called positive) when its deviation exceeds the
    threshold, so TPR/FPR both *increase* as the threshold gets *smaller*
    (stricter detection sensitivity, more things flagged).

    Returns a list of (threshold, tpr, fpr) sorted by fpr ascending.
    """
    pos_scores = [detection_score(r, pos_truth[r["name"]]) for r in pos_scored]
    neg_ok = [r for r in neg_scored if r["status"] == "ok"]
    neg_scores = [detection_score(r, neg_truth[r["name"]]) for r in neg_ok]

    n_pos = len(pos_scored)
    n_neg = len(neg_ok)
    if n_pos == 0 or n_neg == 0:
        return [(0.0, 0.0, 0.0)]

    thresholds = sorted({s for s in pos_scores + neg_scores if s is not None})
    thresholds = [-1.0] + thresholds  # -1: below every real score, so "> t" is true for all of them

    points = []
    for t in thresholds:
        tp = sum(1 for s in pos_scores if s is not None and s > t)
        fp = sum(1 for s in neg_scores if s is not None and s > t)
        points.append((t, tp / n_pos, fp / n_neg))
    points.append((float("inf"), 0.0, 0.0))  # threshold so loose nothing is flagged
    points.sort(key=lambda p: (p[2], p[1]))
    return points


def auc(points):
    """Trapezoidal area under the (fpr, tpr) curve."""
    pts = sorted(set((fpr, tpr) for _, tpr, fpr in points))
    area = 0.0
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        area += (x1 - x0) * (y0 + y1) / 2
    return area


def summarize(scored, group_key=None):
    if group_key:
        groups = defaultdict(list)
        for s in scored:
            groups[s[group_key]].append(s)
    else:
        groups = {"all": scored}

    lines = []
    for key in sorted(groups):
        rows = groups[key]
        n = len(rows)
        n_ok = sum(1 for r in rows if r["status"] == "ok")
        n_pass = sum(1 for r in rows if r["pass"])
        recall = n_pass / n if n else 0.0
        errors = [r["abs_error"] for r in rows if r["status"] == "ok"]
        mean_err = sum(errors) / len(errors) if errors else float("nan")
        lines.append((key, n, n_ok, n_pass, recall, mean_err))
    return lines


def format_table(rows, headers):
    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join(["---"] * len(headers)) + "|"]
    for row in rows:
        lines.append("| " + " | ".join(str(c) for c in row) + " |")
    return "\n".join(lines)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--truth", required=True, help="truth.tsv from select_sv_windows.py")
    p.add_argument("--bedpull-results", nargs="+", required=True, help="one or more bedpull results.tsv (per haplotype)")
    p.add_argument("--liftover-results", nargs="+", required=True, help="one or more liftOver results.tsv (per haplotype)")
    p.add_argument("--paftools-results", nargs="+", help="one or more paftools.js liftover results.tsv (optional extra tool)")
    p.add_argument("--pslmap-results", nargs="+", help="one or more pslMap results.tsv (optional extra tool)")
    p.add_argument("--liftover-multiple-results", nargs="+", help="one or more liftOver -multiple + liftOverMerge results.tsv (optional extra tool)")
    p.add_argument("--bedpull-stitch-results", nargs="+", help="one or more bedpull --stitch_records results.tsv (optional extra tool)")
    p.add_argument("--tolerance-frac", type=float, default=0.05, help="relative length tolerance (default 5%%)")
    p.add_argument("--tolerance-floor", type=float, default=10, help="minimum absolute tolerance in bp (default 10)")
    p.add_argument("--worst-n", type=int, default=20, help="how many worst-offenders to list per tool")
    p.add_argument("--out-report", required=True, help="output markdown report path")
    p.add_argument("--out-details", required=True, help="output per-window TSV path")
    p.add_argument("--neg-truth", help="neg_truth.tsv from select_negative_control_windows.py (enables precision/F1/ROC)")
    p.add_argument("--neg-bedpull-results", nargs="+", help="bedpull results.tsv for the negative-control windows")
    p.add_argument("--neg-liftover-results", nargs="+", help="liftOver results.tsv for the negative-control windows")
    p.add_argument("--neg-paftools-results", nargs="+", help="paftools.js liftover results.tsv for the negative-control windows")
    p.add_argument("--neg-pslmap-results", nargs="+", help="pslMap results.tsv for the negative-control windows")
    p.add_argument("--neg-liftover-multiple-results", nargs="+", help="liftOver -multiple + liftOverMerge results.tsv for the negative-control windows")
    p.add_argument("--neg-bedpull-stitch-results", nargs="+", help="bedpull --stitch_records results.tsv for the negative-control windows")
    p.add_argument("--roc-out", help="output ROC curve plot path; uses kuva if on PATH, else falls back to matplotlib")
    p.add_argument("--roc-data-out", help="output raw (tool, score, label) TSV — kuva/sklearn-ready ROC input")
    args = p.parse_args()

    truth = load_truth(args.truth)

    # (display name, column key, results paths, neg-results paths) — bedpull and
    # liftOver are always present; everything else is an optional extra
    # comparison (see run_paftools_pipeline.sh / run_pslmap_pipeline.sh /
    # run_liftover_multiple_pipeline.sh) that slots in without any other
    # tool-specific code below needing to change.
    tools = [
        {"name": "bedpull", "key": "bedpull", "results": args.bedpull_results, "neg_results": args.neg_bedpull_results},
        {"name": "liftOver", "key": "liftover", "results": args.liftover_results, "neg_results": args.neg_liftover_results},
    ]
    extra_tools = [
        ("paftools.js", "paftools", args.paftools_results, args.neg_paftools_results),
        ("pslMap", "pslmap", args.pslmap_results, args.neg_pslmap_results),
        ("liftOver-multiple", "liftovermultiple", args.liftover_multiple_results, args.neg_liftover_multiple_results),
        ("bedpull+stitch", "bedpullstitch", args.bedpull_stitch_results, args.neg_bedpull_stitch_results),
    ]
    for name, key, results, neg_results in extra_tools:
        if results:
            tools.append({"name": name, "key": key, "results": results, "neg_results": neg_results})

    for tool in tools:
        tool["scored"] = score_tool(truth, load_results(tool["results"]), args.tolerance_frac, args.tolerance_floor)

    have_neg = args.neg_truth and all(tool["neg_results"] for tool in tools)
    if have_neg:
        neg_truth = load_truth(args.neg_truth)
        for tool in tools:
            tool["neg_scored"] = score_tool(neg_truth, load_results(tool["neg_results"]), args.tolerance_frac, args.tolerance_floor)

    # --- per-window details TSV ---
    detail_cols = ["name", "svtype", "svlen", "size_bucket", "expected_length"]
    for tool in tools:
        detail_cols += [f"{tool['key']}_extracted_length", f"{tool['key']}_abs_error", f"{tool['key']}_status", f"{tool['key']}_pass"]
    for tool in tools:
        tool["by_win"] = {r["name"]: r for r in tool["scored"]}
    with open(args.out_details, "w") as f:
        f.write("\t".join(detail_cols) + "\n")
        for name in truth:
            first = tools[0]["by_win"][name]
            values = [name, first["svtype"], first["svlen"], first["size_bucket"], first["expected_length"]]
            for tool in tools:
                r = tool["by_win"][name]
                values += [r["extracted_length"], r["abs_error"], r["status"], r["pass"]]
            f.write("\t".join(str(v) for v in values) + "\n")

    # --- markdown report ---
    title = " vs ".join(tool["name"] for tool in tools) + " — SV extraction accuracy"
    report = []
    report.append(f"# {title}\n")
    report.append(
        f"{len(truth)} SV windows scored. Pass = within "
        f"max({args.tolerance_floor}bp, {args.tolerance_frac:.0%} of expected length).\n"
    )

    roc_curves = {}
    if have_neg:
        report.append(
            f"\n## Precision / Recall / F1\n\n"
            f"Uses {len(neg_truth)} negative-control windows (matched in size, far from any "
            f"known SV — see select_negative_control_windows.py) to give precision a real "
            f"meaning: a false positive is a negative-control window where the tool extracted "
            f"a length that deviates from the no-SV expectation beyond tolerance.\n"
        )
        cm_rows = []
        for tool in tools:
            cm = confusion_matrix(tool["scored"], tool["neg_scored"])
            cm_rows.append(
                (
                    tool["name"],
                    cm["tp"],
                    cm["fp"],
                    cm["fn"],
                    cm["tn"],
                    cm["n_neg_excluded_unmapped"],
                    f"{cm['precision']:.1%}",
                    f"{cm['recall']:.1%}",
                    f"{cm['f1']:.3f}",
                )
            )
            roc_curves[tool["name"]] = roc_points(tool["scored"], truth, tool["neg_scored"], neg_truth)
        report.append(
            format_table(
                cm_rows,
                ["tool", "TP", "FP", "FN", "TN", "neg_excluded(unmapped)", "precision", "recall", "F1"],
            )
        )

        auc_line = "AUC: " + ", ".join(f"{name}={auc(pts):.3f}" for name, pts in roc_curves.items())
        report.append(f"\n{auc_line}\n")

        if args.roc_data_out:
            # Raw (score, label) rows, not pre-aggregated tpr/fpr points — this is
            # what kuva's `roc` subcommand (or sklearn.metrics.roc_curve, or any
            # other ROC tool) actually wants as input; it builds the curve itself.
            with open(args.roc_data_out, "w") as f:
                f.write("tool\tscore\tlabel\n")
                for tool in tools:
                    for r in tool["scored"]:
                        s = detection_score(r, truth[r["name"]])
                        f.write(f"{tool['name']}\t{s if s is not None else -1.0}\t1\n")
                    for r in tool["neg_scored"]:
                        if r["status"] != "ok":
                            continue  # excluded from the confusion matrix too — see module docstring
                        s = detection_score(r, neg_truth[r["name"]])
                        f.write(f"{tool['name']}\t{s}\t0\n")

        if args.roc_out and args.roc_data_out:
            plotted = False
            try:
                subprocess.run(
                    [
                        "kuva",
                        "roc",
                        args.roc_data_out,
                        "--score-col",
                        "score",
                        "--label-col",
                        "label",
                        "--color-by",
                        "tool",
                        "--auc-label",
                        "--title",
                        title.replace("accuracy", "ROC"),
                        "--x-label",
                        "False Positive Rate",
                        "--y-label",
                        "True Positive Rate",
                        "--legend",
                        "Tool",
                        "-o",
                        args.roc_out,
                    ],
                    check=True,
                    capture_output=True,
                    text=True,
                )
                report.append(f"\nROC curve (kuva): `{args.roc_out}`\n")
                plotted = True
            except FileNotFoundError:
                pass  # kuva not on PATH — fall through to the matplotlib fallback
            except subprocess.CalledProcessError as e:
                report.append(f"\n(kuva failed, falling back to matplotlib: {e.stderr.strip()})\n")

            if not plotted:
                try:
                    import matplotlib

                    matplotlib.use("Agg")
                    import matplotlib.pyplot as plt

                    fig, ax = plt.subplots(figsize=(6, 6))
                    colors = {"bedpull": "#2b6cb0", "liftOver": "#c05621", "paftools.js": "#2f855a"}
                    for tool_name, points in roc_curves.items():
                        tprs = [p[1] for p in points]
                        fprs = [p[2] for p in points]
                        ax.plot(
                            fprs,
                            tprs,
                            label=f"{tool_name} (AUC={auc(points):.3f})",
                            color=colors.get(tool_name),
                            linewidth=2,
                        )
                    ax.plot([0, 1], [0, 1], linestyle="--", color="#999999", linewidth=1, label="chance")
                    ax.set_xlim(-0.02, 1.02)
                    ax.set_ylim(-0.02, 1.02)
                    ax.set_xlabel("False Positive Rate")
                    ax.set_ylabel("True Positive Rate")
                    ax.set_title(title.replace("accuracy", "ROC"))
                    ax.legend(loc="lower right")
                    fig.tight_layout()
                    fig.savefig(args.roc_out, dpi=150)
                    plt.close(fig)
                    report.append(f"\nROC curve (matplotlib fallback): `{args.roc_out}`\n")
                except ImportError:
                    report.append(
                        "\n(neither kuva nor matplotlib available — ROC plot skipped; "
                        f"raw scores in `{args.roc_data_out}`)\n"
                    )

    for tool in tools:
        tool_name, scored = tool["name"], tool["scored"]
        report.append(f"\n## {tool_name}\n")
        overall = summarize(scored)
        report.append("### Overall\n")
        report.append(
            format_table(
                [(n, ok, passed, f"{recall:.1%}", f"{err:.1f}") for _, n, ok, passed, recall, err in overall],
                ["n_windows", "n_extracted", "n_pass", "recall", "mean_abs_error(bp)"],
            )
        )
        report.append("\n### By SV type\n")
        by_type = summarize(scored, "svtype")
        report.append(
            format_table(
                [(k, n, ok, passed, f"{recall:.1%}", f"{err:.1f}") for k, n, ok, passed, recall, err in by_type],
                ["svtype", "n_windows", "n_extracted", "n_pass", "recall", "mean_abs_error(bp)"],
            )
        )
        report.append("\n### By size\n")
        by_size = summarize(scored, "size_bucket")
        report.append(
            format_table(
                [(k, n, ok, passed, f"{recall:.1%}", f"{err:.1f}") for k, n, ok, passed, recall, err in by_size],
                ["size_bucket", "n_windows", "n_extracted", "n_pass", "recall", "mean_abs_error(bp)"],
            )
        )
        worst = sorted(
            (r for r in scored if r["status"] == "ok"), key=lambda r: r["abs_error"], reverse=True
        )[: args.worst_n]
        report.append(f"\n### Worst {len(worst)} offenders (extracted but wrong length)\n")
        report.append(
            format_table(
                [
                    (r["name"], r["svtype"], r["svlen"], r["expected_length"], r["extracted_length"], r["abs_error"])
                    for r in worst
                ],
                ["name", "svtype", "svlen", "expected", "extracted", "abs_error"],
            )
        )
        n_unmapped = sum(1 for r in scored if r["status"] == "unmapped")
        report.append(f"\n{n_unmapped} window(s) produced no output at all from {tool_name}.\n")
        n_no_variant = sum(1 for r in scored if r["looks_like_no_variant"])
        if n_no_variant:
            pct = n_no_variant / len(scored)
            report.append(
                f"\nOf the failures, {n_no_variant} ({pct:.1%} of all windows) extracted a length "
                f"matching the *unvaried* window instead — i.e. {tool_name} reported no variant "
                f"here at all, rather than the wrong size. If the truth set and the extraction "
                f"target aren't the same genome (e.g. the hg38↔hs1 confirmation leg), this is "
                f"often the truth set expecting a variant that simply isn't present in the "
                f"sequence being extracted, not a tool failure — see README.md.\n"
            )

    with open(args.out_report, "w") as f:
        f.write("\n".join(report) + "\n")

    print(f"wrote {args.out_report} and {args.out_details}")


if __name__ == "__main__":
    main()
