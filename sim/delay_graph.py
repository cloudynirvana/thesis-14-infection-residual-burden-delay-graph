#!/usr/bin/env python3
"""Toy delayed-risk constraint graph. Research sketch, seed 20260921.

Host-infection and marrow-stress-like evidence are difference constraints.
They are not written into the frozen PK-like vector Theta.
A promotion that would turn waiting times into rate-like coefficients is refused.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SEED = 20260921
ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figures"
OUT = ROOT / "results.json"

# Frozen PK-like symbols. Present so a promotion has a target. Unused by the ranker.
THETA = {"CL": 1.20, "V": 8.00, "ke": 0.15}

# Hypothesis syntax: declared delay bounds. Not estimated. Not in Theta.
# Each edge is (source, target, lo, hi): lo <= t_target - t_source <= hi.
HYPOTHESES = {
    "H0": {
        "name": "burden_chain",
        "nodes": ["S", "B", "R"],
        "edges": [("S", "B", 3.0, 12.0), ("B", "R", 4.0, 20.0)],
    },
    "H1": {
        "name": "infection_serial",
        "nodes": ["S", "I", "B", "R"],
        "edges": [
            ("S", "I", 1.0, 4.0),
            ("I", "B", 3.0, 8.0),
            ("B", "R", 4.0, 20.0),
        ],
    },
    "H2": {
        "name": "competing_host_delays",
        "nodes": ["S", "I", "M", "B", "R"],
        "edges": [
            ("S", "I", 1.0, 4.0),
            ("I", "B", 3.0, 8.0),
            ("B", "M", 2.0, 8.0),
            ("M", "R", 2.0, 10.0),
            ("B", "R", 4.0, 20.0),
        ],
    },
    "H3": {
        "name": "conjunction_before_burden",
        "nodes": ["S", "I", "M", "B", "R"],
        "edges": [
            ("S", "I", 1.0, 4.0),
            ("S", "M", 1.0, 6.0),
            ("I", "B", 2.0, 8.0),
            ("M", "B", 2.0, 8.0),
            ("B", "R", 4.0, 20.0),
        ],
    },
}

# Evidence objects. Store tag is "evidence". The ranker may attach them as
# constraints. It may not copy them into Theta.
SCHEDULES = {
    "burden": [
        {"id": "E_B", "node": "B", "lo": 6.0, "hi": 7.0},
        {"id": "E_R", "node": "R", "lo": 16.0, "hi": 18.0},
    ],
    "host": [
        {"id": "E_I", "node": "I", "lo": 1.5, "hi": 2.5},
        {"id": "E_M", "node": "M", "lo": 10.0, "hi": 12.0},
        {"id": "E_B", "node": "B", "lo": 6.0, "hi": 7.0},
        {"id": "E_R", "node": "R", "lo": 16.0, "hi": 18.0},
    ],
}


def theta_digest(theta: dict) -> str:
    payload = json.dumps(theta, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()


def solve_stn(nodes: list[str], constraints: list[tuple[str, str, float, float]]) -> dict:
    """Simple temporal network. constraints are lo <= t_b - t_a <= hi.

    Node O is the time origin. Shortest paths are the tightest difference bounds.
    A negative cycle means the constraint set is empty.
    """
    names = ["O"] + list(nodes)
    index = {name: i for i, name in enumerate(names)}
    n = len(names)
    dist = [[math.inf] * n for _ in range(n)]
    for i in range(n):
        dist[i][i] = 0.0

    def tighten(src: str, dst: str, weight: float) -> None:
        i, j = index[src], index[dst]
        if weight < dist[i][j]:
            dist[i][j] = weight

    # Pin the origin label S to time 0 when S is in the hypothesis.
    if "S" in index:
        tighten("O", "S", 0.0)
        tighten("S", "O", 0.0)

    for src, dst, lo, hi in constraints:
        tighten(src, dst, hi)
        tighten(dst, src, -lo)

    for k in range(n):
        for i in range(n):
            dik = dist[i][k]
            if dik == math.inf:
                continue
            row_k = dist[k]
            row_i = dist[i]
            for j in range(n):
                cand = dik + row_k[j]
                if cand < row_i[j]:
                    row_i[j] = cand

    for i in range(n):
        if dist[i][i] < -1e-9:
            return {"consistent": False, "windows": {}, "delays": {}}

    windows = {}
    for name in nodes:
        if name == "S":
            windows[name] = [0.0, 0.0]
            continue
        latest = dist[index["O"]][index[name]]
        earliest = -dist[index[name]][index["O"]]
        windows[name] = [round(earliest, 10), round(latest, 10)]

    delays = {}
    for src, dst, lo, hi in constraints:
        if src not in index or dst not in index:
            continue
        tight_hi = dist[index[src]][index[dst]]
        tight_lo = -dist[index[dst]][index[src]]
        key = f"{src}->{dst}"
        delays[key] = [round(tight_lo, 10), round(tight_hi, 10), lo, hi]

    return {"consistent": True, "windows": windows, "delays": delays}


def volume_of(windows: dict) -> float:
    widths = []
    for name, (lo, hi) in windows.items():
        if name == "S":
            continue
        widths.append(max(0.0, hi - lo))
    if not widths:
        return 0.0
    vol = 1.0
    for w in widths:
        vol *= w
    return vol


def evaluate(hyp_id: str, schedule: str) -> dict:
    hyp = HYPOTHESES[hyp_id]
    nodes = hyp["nodes"]
    node_set = set(nodes)
    evidence = SCHEDULES[schedule]
    attached, unexplained = [], []
    for item in evidence:
        if item["node"] in node_set:
            attached.append(item["id"])
        else:
            unexplained.append(item["id"])
    constraints = list(hyp["edges"])
    for item in evidence:
        if item["node"] in node_set:
            constraints.append(("O", item["node"], item["lo"], item["hi"]))
    solved = solve_stn(nodes, constraints)
    vol = volume_of(solved["windows"]) if solved["consistent"] else None
    return {
        "id": hyp_id,
        "name": hyp["name"],
        "schedule": schedule,
        "consistent": solved["consistent"],
        "unexplained": unexplained,
        "n_unexplained": len(unexplained),
        "attached": attached,
        "n_nodes": len(nodes),
        "volume": None if vol is None else round(vol, 10),
        "windows": solved["windows"],
        "delays": solved["delays"],
    }


def lex_key(row: dict) -> tuple:
    # Pre-declared. Unhosted evidence outranks a narrow window.
    if not row["consistent"]:
        return (1, row["n_unexplained"], math.inf, row["n_nodes"], row["id"])
    return (0, row["n_unexplained"], row["volume"], row["n_nodes"], row["id"])


def volume_key(row: dict) -> tuple:
    # Sensitivity: ignore unexplained evidence. Not the ranker used in Chapter Four.
    if not row["consistent"]:
        return (1, math.inf, row["id"])
    return (0, row["volume"], row["id"])


def rank_schedule(schedule: str) -> dict:
    rows = [evaluate(hid, schedule) for hid in HYPOTHESES]
    ordered = sorted(rows, key=lex_key)
    by_volume = sorted(rows, key=volume_key)
    for i, row in enumerate(ordered, start=1):
        row["rank"] = i
    return {
        "schedule": schedule,
        "order": [row["id"] for row in ordered],
        "volume_only_order": [row["id"] for row in by_volume],
        "rows": ordered,
    }


def propose_pk_coefficients(theta: dict) -> dict:
    """Illegal map from evidence midpoints onto rate-like symbols.

    The numbers are computed so the refusal transcript can show them.
    They are not estimates and are not written back.
    """
    by_id = {item["id"]: item for item in SCHEDULES["host"]}
    t_i = 0.5 * (by_id["E_I"]["lo"] + by_id["E_I"]["hi"])
    t_m = 0.5 * (by_id["E_M"]["lo"] + by_id["E_M"]["hi"])
    t_b = 0.5 * (by_id["E_B"]["lo"] + by_id["E_B"]["hi"])
    t_r = 0.5 * (by_id["E_R"]["lo"] + by_id["E_R"]["hi"])
    wait_host = t_m - t_b
    wait_relapse = t_r - t_b
    k_inf = 1.0 / t_i
    k_host = 1.0 / wait_host
    # V is read and not written. CL and ke would be overwritten. k_inf and k_host are new keys.
    return {
        "source_store": "evidence",
        "source_ids": ["E_I", "E_M", "E_B", "E_R"],
        "attempted_write": {
            "k_inf": k_inf,
            "k_host": k_host,
            "CL": k_host * theta["V"],
            "ke": math.log(2.0) / wait_relapse,
        },
        "midpoints": {"t_I": t_i, "t_M": t_m, "t_B": t_b, "t_R": t_r},
    }


def refuse_promotion(theta: dict) -> dict:
    """Refuse any write whose source store is evidence or whose keys leave Theta."""
    proposal = propose_pk_coefficients(theta)
    allowed = set(THETA)
    new_keys = sorted(set(proposal["attempted_write"]) - allowed)
    overwritten = sorted(set(proposal["attempted_write"]) & allowed)
    reasons = []
    if proposal["source_store"] == "evidence":
        reasons.append("source store is evidence, which is a non-parameter")
    if new_keys:
        reasons.append("new keys are PK-like coefficients: " + ", ".join(new_keys))
    if overwritten:
        reasons.append("existing Theta keys would be overwritten from evidence: " + ", ".join(overwritten))
    # The ranker is not called. Theta is returned as the same object contents.
    theta_after = dict(theta)
    return {
        "status": "REFUSED",
        "code": "ILLEGAL_PROMOTION",
        "reasons": reasons,
        "proposal_not_a_result": proposal,
        "theta_before": dict(theta),
        "theta_after": theta_after,
        "digest_before": theta_digest(theta),
        "digest_after": theta_digest(theta_after),
        "ranker_called": False,
    }


def sample_windows(row: dict, n: int, rng: np.random.Generator) -> dict:
    """Rejection sample inside propagated windows, then enforce delay edges."""
    if not row["consistent"]:
        return {"accepted": 0, "tried": 0, "points": []}
    hyp = HYPOTHESES[row["id"]]
    windows = row["windows"]
    free = [name for name in hyp["nodes"] if name != "S"]
    accepted = []
    tried = 0
    limit = n * 400
    while len(accepted) < n and tried < limit:
        tried += 1
        point = {"S": 0.0}
        ok = True
        for name in free:
            lo, hi = windows[name]
            if hi < lo:
                ok = False
                break
            point[name] = float(rng.uniform(lo, hi)) if hi > lo else lo
        if not ok:
            break
        for src, dst, lo, hi in hyp["edges"]:
            delay = point[dst] - point[src]
            if delay < lo - 1e-9 or delay > hi + 1e-9:
                ok = False
                break
        if ok:
            accepted.append({k: round(point[k], 6) for k in hyp["nodes"]})
    return {"accepted": len(accepted), "tried": tried, "points": accepted}


def plot_topologies(path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.6))
    specs = [
        (axes[0], "H0", ["S", "B", "R"], [("S", "B"), ("B", "R")]),
        (
            axes[1],
            "H2",
            ["S", "I", "B", "M", "R"],
            [("S", "I"), ("I", "B"), ("B", "M"), ("M", "R"), ("B", "R")],
        ),
    ]
    pos = {
        "S": (0.0, 0.0),
        "I": (1.3, 0.85),
        "B": (2.6, 0.0),
        "M": (3.9, 0.85),
        "R": (5.2, 0.0),
    }
    for ax, hid, nodes, edges in specs:
        hyp = HYPOTHESES[hid]
        bounds = {(a, b): (lo, hi) for a, b, lo, hi in hyp["edges"]}
        for src, dst in edges:
            x0, y0 = pos[src]
            x1, y1 = pos[dst]
            ax.annotate(
                "",
                xy=(x1, y1),
                xytext=(x0, y0),
                arrowprops=dict(arrowstyle="-|>", color="#1f4e79", lw=1.2),
            )
            lo, hi = bounds[(src, dst)]
            label = f"[{lo:g}, {hi:g}]"
            ax.text((x0 + x1) / 2, (y0 + y1) / 2 + 0.16, label, ha="center", va="bottom", fontsize=8, color="#333333")
        for name in nodes:
            x, y = pos[name]
            ax.scatter([x], [y], s=280, c="#f7f1e8", edgecolors="#1f4e79", zorder=3)
            ax.text(x, y, name, ha="center", va="center", fontsize=10, zorder=4)
        ax.set_title(f"{hid}  {hyp['name'].replace('_', ' ')}", fontsize=10)
        ax.set_xlim(-0.6, 5.8)
        ax.set_ylim(-0.7, 1.55)
        ax.axis("off")
    fig.suptitle("Declared delay bounds, not fitted rates", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_windows(burden_rows: dict, host_rows: dict, path: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(8.4, 5.4), sharex=True)
    panels = [
        (axes[0], "Burden evidence only", burden_rows),
        (axes[1], "Host-infection evidence attached", host_rows),
    ]
    colors = {"H0": "#1f4e79", "H2": "#b85c38"}
    for ax, title, rows in panels:
        ypos = 0
        yticks, ylabels = [], []
        for hid in ("H0", "H2"):
            row = rows[hid]
            for name in ["I", "B", "M", "R"]:
                yticks.append(ypos)
                ylabels.append(f"{hid}  {name}")
                if name not in row["windows"]:
                    ax.plot([0], [ypos], marker="x", color="#999999")
                    ax.text(0.15, ypos, "no node", va="center", fontsize=8, color="#666666")
                else:
                    lo, hi = row["windows"][name]
                    ax.hlines(ypos, lo, hi, colors=colors[hid], lw=3)
                    ax.plot([lo, hi], [ypos, ypos], "o", color=colors[hid], ms=4)
                ypos += 1
            ypos += 0.4
        ax.set_yticks(yticks)
        ax.set_yticklabels(ylabels, fontsize=8)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("Model time")
        ax.grid(axis="x", color="#dddddd")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_ranks(burden: dict, host: dict, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    order_a = {hid: i for i, hid in enumerate(burden["order"])}
    order_b = {hid: i for i, hid in enumerate(host["order"])}
    colors = {"H0": "#1f4e79", "H1": "#3d6b4f", "H2": "#b85c38", "H3": "#6b6b6b"}
    for hid in HYPOTHESES:
        ax.plot([0, 1], [order_a[hid], order_b[hid]], color=colors[hid], lw=2)
        ax.scatter([0, 1], [order_a[hid], order_b[hid]], color=colors[hid], s=40, zorder=3)
        ax.text(-0.06, order_a[hid], hid, ha="right", va="center", color=colors[hid], fontsize=10)
        ax.text(1.06, order_b[hid], hid, ha="left", va="center", color=colors[hid], fontsize=10)
    ax.set_xlim(-0.35, 1.35)
    ax.set_ylim(3.4, -0.4)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Burden evidence", "Host evidence attached"])
    ax.set_ylabel("Rank (1 is preferred)")
    ax.set_yticks([0, 1, 2, 3])
    ax.set_yticklabels(["1", "2", "3", "4"])
    ax.set_title("Lexicographic rank. Theta is not an input.")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_refusal(refusal: dict, path: Path) -> None:
    before = refusal["theta_before"]
    attempted = refusal["proposal_not_a_result"]["attempted_write"]
    rows = [
        ["CL", f"{before['CL']:.2f}", f"{attempted['CL']:.4f}  refused"],
        ["V", f"{before['V']:.2f}", "not written"],
        ["ke", f"{before['ke']:.2f}", f"{attempted['ke']:.5f}  refused"],
        ["k_inf", "absent from Theta", f"{attempted['k_inf']:.4f}  refused"],
        ["k_host", "absent from Theta", f"{attempted['k_host']:.4f}  refused"],
    ]
    fig, ax = plt.subplots(figsize=(7.8, 3.2))
    ax.axis("off")
    ax.set_title("ILLEGAL_PROMOTION. The right-hand column is not a result.", fontsize=11, pad=8)
    table = ax.table(
        cellText=rows,
        colLabels=["Symbol", "Theta", "Proposal"],
        loc="center",
        cellLoc="left",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.0, 1.45)
    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor("#555555")
        if r == 0:
            cell.set_facecolor("#eeeeee")
            cell.set_text_props(weight="bold")
        elif c == 2 and r != 0:
            cell.set_facecolor("#f4f4f4")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    burden = rank_schedule("burden")
    host = rank_schedule("host")
    refusal = refuse_promotion(THETA)

    burden_by = {row["id"]: row for row in burden["rows"]}
    host_by = {row["id"]: row for row in host["rows"]}
    samples = {
        "burden_H0": sample_windows(burden_by["H0"], 200, rng),
        "host_H2": sample_windows(host_by["H2"], 200, rng),
    }
    # Points are for the figure audit trail; keep a short prefix in JSON.
    sample_summary = {
        key: {"accepted": val["accepted"], "tried": val["tried"], "head": val["points"][:5]}
        for key, val in samples.items()
    }

    plot_topologies(FIG / "topologies.png")
    plot_windows(burden_by, host_by, FIG / "windows.png")
    plot_ranks(burden, host, FIG / "rank_change.png")
    plot_refusal(refusal, FIG / "refusal.png")

    # Drop bulky point clouds from the stored rows (windows and delays stay).
    def slim(block: dict) -> dict:
        return {
            "schedule": block["schedule"],
            "order": block["order"],
            "volume_only_order": block["volume_only_order"],
            "rows": block["rows"],
        }

    results = {
        "seed": SEED,
        "theta": THETA,
        "theta_digest": theta_digest(THETA),
        "hypotheses": {
            hid: {"name": hyp["name"], "nodes": hyp["nodes"], "edges": hyp["edges"]}
            for hid, hyp in HYPOTHESES.items()
        },
        "schedules": SCHEDULES,
        "rank_key": [
            "infeasible last",
            "fewer unexplained evidence objects",
            "smaller feasible-window volume",
            "fewer nodes",
            "hypothesis id",
        ],
        "burden": slim(burden),
        "host": slim(host),
        "refusal": refusal,
        "samples": sample_summary,
        "notes": [
            "Evidence windows are declared. They are not a cohort.",
            "Structural delay bounds are hypothesis syntax, not Theta.",
            "The refused coefficients are not estimates.",
        ],
    }
    OUT.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")

    assert refusal["status"] == "REFUSED"
    assert refusal["digest_before"] == refusal["digest_after"]
    assert refusal["theta_before"] == refusal["theta_after"] == THETA
    assert refusal["ranker_called"] is False
    assert burden["order"][0] == "H0"
    assert host["order"][0] == "H2"
    assert burden["order"] != host["order"]
    assert host_by["H3"]["consistent"] is False
    assert host_by["H0"]["n_unexplained"] == 2
    assert host_by["H2"]["n_unexplained"] == 0
    assert "k_inf" not in THETA and "k_host" not in THETA
    print(json.dumps({
        "burden_order": burden["order"],
        "host_order": host["order"],
        "burden_volume_only": burden["volume_only_order"],
        "host_volume_only": host["volume_only_order"],
        "refusal": refusal["code"],
        "digest": refusal["digest_before"],
    }, indent=2))


if __name__ == "__main__":
    main()
