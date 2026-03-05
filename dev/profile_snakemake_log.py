#!/usr/bin/env python3
"""Parse a Snakemake log file and produce an execution profile report.

Usage:
    python dev/profile_snakemake_log.py .snakemake/log/<logfile>.snakemake.log
    python dev/profile_snakemake_log.py --html report.html <logfile>
    python dev/profile_snakemake_log.py --latest          # auto-pick newest log
"""

import argparse
import re
import sys
from collections import defaultdict
from datetime import datetime, timedelta
from pathlib import Path

# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

TIMESTAMP_RE = re.compile(r"^\[(.+?)\]$")
RULE_START_RE = re.compile(r"^(?:localrule|rule) (\S+):$")
JOBID_RE = re.compile(r"^\s+jobid:\s*(\d+)")
WILDCARDS_RE = re.compile(r"^\s+wildcards:\s*(.+)")
THREADS_RE = re.compile(r"^\s+threads:\s*(\d+)")
FINISHED_RE = re.compile(r"^Finished jobid:\s*(\d+)\s*\(Rule:\s*(\S+)\)")
STEPS_RE = re.compile(r"^(\d+) of (\d+) steps")


def parse_log(path):
    """Parse Snakemake log, return list of completed job records."""
    jobs_pending = {}  # jobid -> {rule, start, wildcards, threads}
    jobs_done = []
    current_ts = None
    current_rule = None

    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")

            m = TIMESTAMP_RE.match(line)
            if m:
                try:
                    current_ts = datetime.strptime(m.group(1), "%a %b %d %H:%M:%S %Y")
                except ValueError:
                    pass
                continue

            m = RULE_START_RE.match(line)
            if m:
                current_rule = m.group(1)
                continue

            m = JOBID_RE.match(line)
            if m and current_rule and current_ts:
                jobid = int(m.group(1))
                jobs_pending[jobid] = {
                    "rule": current_rule,
                    "start": current_ts,
                    "wildcards": "",
                    "threads": 1,
                }
                current_rule = None
                continue

            m = WILDCARDS_RE.match(line)
            if m:
                # attach to most recently seen jobid
                for jid in reversed(list(jobs_pending)):
                    if "wildcards" in jobs_pending[jid]:
                        jobs_pending[jid]["wildcards"] = m.group(1)
                        break
                continue

            m = THREADS_RE.match(line)
            if m:
                for jid in reversed(list(jobs_pending)):
                    jobs_pending[jid]["threads"] = int(m.group(1))
                    break
                continue

            m = FINISHED_RE.match(line)
            if m and current_ts:
                jobid = int(m.group(1))
                if jobid in jobs_pending:
                    rec = jobs_pending.pop(jobid)
                    rec["end"] = current_ts
                    rec["duration"] = (rec["end"] - rec["start"]).total_seconds()
                    rec["jobid"] = jobid
                    jobs_done.append(rec)
                continue

    return jobs_done


# ---------------------------------------------------------------------------
# Analysis helpers
# ---------------------------------------------------------------------------

def fmt_duration(seconds):
    """Format seconds as human-readable duration."""
    if seconds < 60:
        return f"{seconds:.0f}s"
    elif seconds < 3600:
        m, s = divmod(seconds, 60)
        return f"{int(m)}m {int(s)}s"
    else:
        h, rem = divmod(seconds, 3600)
        m, s = divmod(rem, 60)
        return f"{int(h)}h {int(m)}m {int(s)}s"


def analyze(jobs):
    """Compute summary statistics from job records."""
    if not jobs:
        return {}

    total_start = min(j["start"] for j in jobs)
    total_end = max(j["end"] for j in jobs)
    wall_time = (total_end - total_start).total_seconds()

    # Per-rule aggregates
    by_rule = defaultdict(lambda: {"count": 0, "total_sec": 0.0, "max_sec": 0.0, "jobs": []})
    for j in jobs:
        r = by_rule[j["rule"]]
        r["count"] += 1
        r["total_sec"] += j["duration"]
        r["max_sec"] = max(r["max_sec"], j["duration"])
        r["jobs"].append(j)

    # Sort by total time descending
    rules_ranked = sorted(by_rule.items(), key=lambda x: -x[1]["total_sec"])

    # Top 10 longest individual jobs
    top_jobs = sorted(jobs, key=lambda j: -j["duration"])[:15]

    # Phase analysis: group rules into logical phases
    phase_map = {
        "get_fastq": "Download", "process_fastq": "Download",
        "filter_chip": "Alignment", "filter_rna": "Alignment",
        "STAR_map": "Alignment", "map_chromap": "Alignment",
        "bowtie2_map": "Alignment",
        "filter_structural_rna": "Alignment", "dispatch_srna_fastq": "Alignment",
        "merge_replicates": "Merge", "making_pseudo_replicates": "Merge",
        "calling_peaks": "Peak calling", "best_peaks": "Peak calling",
        "idr_analysis": "IDR", "select_peaks": "IDR",
        "perform_pairwise_diff": "Differential",
        "call_all_DEGs": "Differential", "call_all_differential": "Differential",
        "make_bigwig": "Tracks", "make_rna_stranded": "Tracks",
        "make_srna_stranded": "Tracks",
        "computing_matrix": "Heatmaps", "plot_heatmap": "Heatmaps",
        "combined_analysis": "Heatmaps",
        "make_fingerprint": "QC", "make_chip_stats": "QC",
        "make_rna_stats": "QC", "make_peak_stats": "QC",
        "make_srna_size": "QC", "has_header": "QC", "is_stranded": "QC",
        "filter_size_srna": "sRNA analysis",
        "analyze_all_srna": "sRNA analysis",
        "make_cluster": "sRNA analysis", "combine_cluster": "sRNA analysis",
    }

    phases = defaultdict(lambda: {"total_sec": 0.0, "count": 0})
    for j in jobs:
        phase = "Other"
        for prefix, p in phase_map.items():
            if j["rule"].startswith(prefix):
                phase = p
                break
        phases[phase]["total_sec"] += j["duration"]
        phases[phase]["count"] += j["count"] if "count" in j else 1

    return {
        "total_jobs": len(jobs),
        "wall_time": wall_time,
        "total_start": total_start,
        "total_end": total_end,
        "cpu_time": sum(j["duration"] * j["threads"] for j in jobs),
        "rules_ranked": rules_ranked,
        "top_jobs": top_jobs,
        "phases": dict(sorted(phases.items(), key=lambda x: -x[1]["total_sec"])),
    }


# ---------------------------------------------------------------------------
# Markdown report
# ---------------------------------------------------------------------------

def report_markdown(stats):
    lines = []
    lines.append("# Snakemake Execution Profile\n")
    lines.append(f"**Run period:** {stats['total_start'].strftime('%Y-%m-%d %H:%M:%S')} — "
                 f"{stats['total_end'].strftime('%H:%M:%S')}")
    lines.append(f"**Wall time:** {fmt_duration(stats['wall_time'])}")
    lines.append(f"**Total jobs:** {stats['total_jobs']}")
    lines.append(f"**Total CPU time:** {fmt_duration(stats['cpu_time'])}")
    lines.append(f"**Parallelism:** {stats['cpu_time'] / stats['wall_time']:.1f}x average\n")

    # Phase summary
    lines.append("## Phase Summary\n")
    lines.append("| Phase | Jobs | Total time | % of wall |")
    lines.append("|-------|------|-----------|-----------|")
    for phase, data in stats["phases"].items():
        pct = data["total_sec"] / stats["wall_time"] * 100
        lines.append(f"| {phase} | {data['count']} | "
                     f"{fmt_duration(data['total_sec'])} | {pct:.1f}% |")

    # Per-rule table
    lines.append("\n## Per-Rule Breakdown\n")
    lines.append("| Rule | Count | Total | Mean | Max |")
    lines.append("|------|-------|-------|------|-----|")
    for rule, data in stats["rules_ranked"]:
        mean = data["total_sec"] / data["count"]
        lines.append(f"| {rule} | {data['count']} | "
                     f"{fmt_duration(data['total_sec'])} | "
                     f"{fmt_duration(mean)} | "
                     f"{fmt_duration(data['max_sec'])} |")

    # Top individual jobs
    lines.append("\n## Slowest Individual Jobs\n")
    lines.append("| # | Rule | Duration | Wildcards |")
    lines.append("|---|------|----------|-----------|")
    for i, j in enumerate(stats["top_jobs"], 1):
        wc = j["wildcards"][:60] if j["wildcards"] else ""
        lines.append(f"| {i} | {j['rule']} | {fmt_duration(j['duration'])} | {wc} |")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# HTML report
# ---------------------------------------------------------------------------

def report_html(stats):
    md = report_markdown(stats)
    # Convert simple markdown tables and headers to HTML
    rows = md.split("\n")
    body_parts = []
    in_table = False

    for row in rows:
        if row.startswith("# "):
            if in_table:
                body_parts.append("</table>")
                in_table = False
            body_parts.append(f"<h1>{row[2:]}</h1>")
        elif row.startswith("## "):
            if in_table:
                body_parts.append("</table>")
                in_table = False
            body_parts.append(f"<h2>{row[3:]}</h2>")
        elif row.startswith("**"):
            body_parts.append(f"<p>{row}</p>")
        elif row.startswith("|") and not row.startswith("|--"):
            cells = [c.strip() for c in row.split("|")[1:-1]]
            if not in_table:
                body_parts.append('<table>')
                tag = "th"
                in_table = True
            else:
                tag = "td"
            body_parts.append(
                "<tr>" + "".join(f"<{tag}>{c}</{tag}>" for c in cells) + "</tr>"
            )
        elif row.startswith("|--"):
            continue  # skip markdown separator
        elif row.strip() == "":
            if in_table:
                body_parts.append("</table>")
                in_table = False
        else:
            body_parts.append(f"<p>{row}</p>")

    if in_table:
        body_parts.append("</table>")

    # Gantt-style timeline chart (SVG)
    gantt = _build_gantt_svg(stats)

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Snakemake Execution Profile</title>
<style>
  body {{ font-family: -apple-system, system-ui, sans-serif; max-width: 1200px;
         margin: 2rem auto; padding: 0 1rem; color: #1a1a1a; }}
  h1 {{ border-bottom: 2px solid #2563eb; padding-bottom: 0.3rem; }}
  h2 {{ color: #2563eb; margin-top: 2rem; }}
  table {{ border-collapse: collapse; width: 100%; margin: 0.5rem 0 1.5rem; }}
  th {{ background: #2563eb; color: white; text-align: left; padding: 0.5rem 0.75rem; }}
  td {{ padding: 0.4rem 0.75rem; border-bottom: 1px solid #e5e7eb; }}
  tr:hover td {{ background: #f0f4ff; }}
  p {{ line-height: 1.6; }}
  .gantt {{ margin: 1rem 0; overflow-x: auto; }}
  svg text {{ font-size: 11px; font-family: monospace; }}
</style>
</head>
<body>
{"".join(body_parts)}
<h2>Execution Timeline</h2>
<div class="gantt">{gantt}</div>
</body>
</html>"""


def _build_gantt_svg(stats):
    """Build a simple SVG Gantt chart of job execution."""
    jobs = []
    for rule, data in stats["rules_ranked"]:
        for j in data["jobs"]:
            jobs.append(j)

    if not jobs:
        return "<p>No job data available.</p>"

    t0 = stats["total_start"]
    wall = stats["wall_time"]
    if wall == 0:
        return "<p>Zero wall time.</p>"

    # Assign colors by rule
    rules = list(dict.fromkeys(j["rule"] for j in sorted(jobs, key=lambda x: x["start"])))
    palette = [
        "#2563eb", "#dc2626", "#16a34a", "#9333ea", "#ea580c",
        "#0891b2", "#c026d3", "#65a30d", "#d97706", "#475569",
        "#e11d48", "#059669", "#7c3aed", "#db2777", "#0d9488",
    ]
    rule_color = {r: palette[i % len(palette)] for i, r in enumerate(rules)}

    # Layout: rows for concurrent jobs
    row_height = 18
    label_width = 200
    chart_width = 700
    svg_width = label_width + chart_width + 20

    # Assign rows: greedy lane assignment
    lanes = []  # each lane has an end_time
    job_rows = []
    sorted_jobs = sorted(jobs, key=lambda j: j["start"])
    for j in sorted_jobs:
        s = (j["start"] - t0).total_seconds()
        placed = False
        for li, lane_end in enumerate(lanes):
            if s >= lane_end:
                lanes[li] = (j["end"] - t0).total_seconds()
                job_rows.append((j, li))
                placed = True
                break
        if not placed:
            lanes.append((j["end"] - t0).total_seconds())
            job_rows.append((j, len(lanes) - 1))

    num_lanes = len(lanes)
    svg_height = max(num_lanes * row_height + 40, 80)

    parts = [f'<svg width="{svg_width}" height="{svg_height}" xmlns="http://www.w3.org/2000/svg">']

    # Time axis ticks
    for pct in range(0, 101, 25):
        x = label_width + (pct / 100) * chart_width
        t_sec = pct / 100 * wall
        parts.append(f'<line x1="{x}" y1="0" x2="{x}" y2="{svg_height - 20}" '
                     f'stroke="#e5e7eb" stroke-width="1"/>')
        parts.append(f'<text x="{x}" y="{svg_height - 5}" text-anchor="middle" '
                     f'fill="#6b7280">{fmt_duration(t_sec)}</text>')

    # Job bars
    for j, lane in job_rows:
        s = (j["start"] - t0).total_seconds()
        e = (j["end"] - t0).total_seconds()
        x = label_width + (s / wall) * chart_width
        w = max(((e - s) / wall) * chart_width, 2)
        y = lane * row_height + 2
        h = row_height - 4
        color = rule_color.get(j["rule"], "#475569")
        wc = j["wildcards"][:40] if j["wildcards"] else ""
        title = f'{j["rule"]} ({fmt_duration(j["duration"])})\n{wc}'
        parts.append(f'<rect x="{x:.1f}" y="{y}" width="{w:.1f}" height="{h}" '
                     f'fill="{color}" rx="2" opacity="0.85">'
                     f'<title>{title}</title></rect>')

    # Legend (top rules)
    legend_y = svg_height - 18
    # Put legend below, extend SVG
    legend_height = ((len(rules) - 1) // 5 + 1) * 18 + 10
    parts[0] = parts[0].replace(f'height="{svg_height}"',
                                f'height="{svg_height + legend_height}"')

    for i, rule in enumerate(rules[:15]):
        lx = (i % 5) * 230 + 10
        ly = svg_height + (i // 5) * 18 + 12
        color = rule_color[rule]
        parts.append(f'<rect x="{lx}" y="{ly - 10}" width="12" height="12" fill="{color}" rx="2"/>')
        parts.append(f'<text x="{lx + 16}" y="{ly}" fill="#1a1a1a">{rule}</text>')

    parts.append("</svg>")
    return "\n".join(parts)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def find_latest_log():
    log_dir = Path(".snakemake/log")
    if not log_dir.exists():
        sys.exit("No .snakemake/log directory found. Run from the Snakemake project root.")
    logs = sorted(log_dir.glob("*.snakemake.log"))
    if not logs:
        sys.exit("No Snakemake log files found.")
    return logs[-1]


def main():
    parser = argparse.ArgumentParser(description="Profile a Snakemake run from its log file.")
    parser.add_argument("logfile", nargs="?", help="Path to .snakemake.log file")
    parser.add_argument("--latest", action="store_true", help="Auto-select newest log")
    parser.add_argument("--html", metavar="FILE", help="Write HTML report to FILE")
    args = parser.parse_args()

    if args.latest or not args.logfile:
        logpath = find_latest_log()
    else:
        logpath = Path(args.logfile)

    if not logpath.exists():
        sys.exit(f"Log file not found: {logpath}")

    print(f"Parsing {logpath} ...", file=sys.stderr)
    jobs = parse_log(logpath)

    if not jobs:
        sys.exit("No completed jobs found in log.")

    stats = analyze(jobs)

    if args.html:
        html = report_html(stats)
        Path(args.html).write_text(html)
        print(f"HTML report written to {args.html}", file=sys.stderr)
    else:
        print(report_markdown(stats))


if __name__ == "__main__":
    main()
