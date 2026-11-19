"""Finite-window occupancy bounds for a two-threshold piecewise system."""
import argparse
from hashlib import sha256
from itertools import combinations
import json
from math import log, tanh
from pathlib import Path
from random import Random
from statistics import mean, median

PROTOCOL_HASH = "939129ff41d9d6fa82703af6b39e3a13ff999d0ee59b4e8eff3c6151a97aa9eb"
HORIZON = 32
PUBLISHED = (
    ((2, 1, 3, 3), (2, 2, 1, 4)),
    ((1, 2, 2, 1), (3, 2, 3, 4)),
    ((3, 2, 3, 2), (1, 1, 4, 4)),
)
KERNELS = ("instant_sign", "hysteresis", "tanh")
METRICS = ((0, 0, 0, 1), (0, 1, 0, 1), (0, 0, 1, 1))


def coefficients(a, b):
    return (a[0], -a[1], a[2], -a[3]), (b[0], b[1], -b[2], -b[3])


def trial(a, b, initial, dt, width, kernel):
    vx, vy = coefficients(a, b)
    x, y = initial
    u, v = int(x >= 0), int(y >= 0)
    counts = [0.0] * 4
    steps = round(HORIZON / dt)
    for _ in range(steps):
        if kernel == "tanh":
            u = (1 + tanh(x / width)) / 2
            v = (1 + tanh(y / width)) / 2
            p = ((1-u)*(1-v), u*(1-v), (1-u)*v, u*v)
        else:
            if kernel == "instant_sign":
                u, v = int(x >= 0), int(y >= 0)
            else:
                if abs(x) >= width:
                    u = int(x >= 0)
                if abs(y) >= width:
                    v = int(y >= 0)
            index = int(u) + 2 * int(v)
            p = tuple(int(i == index) for i in range(4))
        for i in range(4):
            counts[i] += p[i]
        x += dt * dot(p, vx)
        y += dt * dot(p, vy)
    return tuple(c / steps for c in counts), (x, y)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def solved(rows, targets):
    """Solve a four-variable linear system. Return None if rank deficient."""
    table = [list(row) + [t] for row, t in zip(rows, targets)]
    for i in range(4):
        pivot = max(range(i, 4), key=lambda j: abs(table[j][i]))
        if abs(table[pivot][i]) < 1e-12:
            return None
        table[i], table[pivot] = table[pivot], table[i]
        denominator = table[i][i]
        table[i] = [x / denominator for x in table[i]]
        for j in range(4):
            if j != i:
                multiple = table[j][i]
                table[j] = [y - multiple * x for x, y in zip(table[i], table[j])]
    return tuple(row[-1] for row in table)


def occupancy_vertices(vx, vy, x_band, y_band):
    inequalities = [(tuple(-int(i == j) for i in range(4)), 0) for j in range(4)]
    for vector, (lower, upper) in [(vx, x_band), (vy, y_band)]:
        inequalities.append((vector, upper))
        inequalities.append((tuple(-i for i in vector), -lower))
    vertices = []
    for selected in combinations(inequalities, 3):
        solution = solved([(1, 1, 1, 1), *(r for r, _ in selected)],
                          [1, *(n for _, n in selected)])
        if solution is not None and all(dot(c, solution) <= d + 1e-9
                                        for c, d in inequalities):
            vertices.append(solution)
    assert vertices
    return vertices


def bounds(a, b, state, dt, width):
    vx, vy = coefficients(a, b)
    # Once beyond each threshold, every admissible rule pushes inward.
    def cap(current, rates):
        threshold = max(width, 0.5 * width * log(max(rates) / min(rates)))
        return max(abs(current), threshold + dt * max(rates))
    xmax, ymax = (cap(x, r) for x, r in zip(state, (a, b)))
    vertices = occupancy_vertices(vx, vy, ((-xmax-state[0])/HORIZON, (xmax-state[0])/HORIZON),
                                  ((-ymax-state[1])/HORIZON, (ymax-state[1])/HORIZON))
    intervals = []
    for metric in METRICS:
        vals = [dot(metric, vertex) for vertex in vertices]
        intervals.append((max(0.0, min(vals)), min(1.0, max(vals))))
    return intervals, (xmax, ymax)


def maximum_entropy(vx, vy):
    vertices = occupancy_vertices(vx, vy, (0, 0), (0, 0))
    start, end = min(vertices), max(vertices)
    blend = lambda z: tuple((1-z)*a+z*b for a,b in zip(start,end))
    def entropy(z):
        return -sum(p*log(p) for p in blend(z) if p > 0)
    lo, hi = 0.0, 1.0
    for _ in range(70):
        p, q = (2*lo + hi)/3, (lo + 2*hi)/3
        if entropy(p) < entropy(q):
            lo = p
        else:
            hi = q
    return blend((lo+hi)/2)


def product_closure(a, b):
    def p_of(q):
        inc, dec = (1-q)*a[0]+q*a[2], (1-q)*a[1]+q*a[3]
        return inc / (inc+dec)
    def q_of(p):
        inc, dec = (1-p)*b[0]+p*b[1], (1-p)*b[2]+p*b[3]
        return inc / (inc+dec)
    lo, hi = 0.0, 1.0
    for _ in range(64):
        q = (lo+hi)/2
        if q < q_of(p_of(q)):
            lo = q
        else:
            hi = q
    q = (lo+hi)/2
    p = p_of(q)
    return ((1-p)*(1-q), p*(1-q), (1-p)*q, p*q)


def models():
    rng = Random(6112026)
    matched = [(tuple(rng.randint(1, 5) for _ in range(4)),
                tuple(rng.randint(1, 5) for _ in range(4))) for _ in range(24)]
    affine = []
    for _ in range(16):
        a, b, c, d = (rng.randint(1, 5) for _ in range(4))
        affine.append(((a,b,a,b),(c,c,d,d)))
    return {"published_coefficients": PUBLISHED, "generated_coupled": matched,
            "generated_uncoupled": affine}


def evaluate(field, key):
    a, b = field
    vx, vy = coefficients(a, b)
    assert abs(dot(vx, product_closure(a,b))) < 1e-10
    assert abs(dot(vy, product_closure(a,b))) < 1e-10
    product = product_closure(a, b)
    entropy = maximum_entropy(vx, vy)
    rng = Random(70200 + key)
    initial = [tuple(rng.uniform(-.18,.18) for _ in range(2)) for _ in range(2)]
    errors = {"product_11": [], "entropy_11": [], "unconstrained_interval_11": [],
              "balance_11_width": [], "balance_first_width": [], "balance_second_width": []}
    ranges, max_telescope, covered, trials = [], 0.0, [0, 0, 0], 0
    for state in initial:
        for dt in (.10, .04):
            for multiplier in (.5, 2.0):
                width = multiplier * dt
                intervals, caps = bounds(a, b, state, dt, width)
                cohort = []
                for kernel in KERNELS:
                    occupancy, terminal = trial(a, b, state, dt, width, kernel)
                    assert abs(sum(occupancy)-1) < 1e-10
                    max_telescope = max(max_telescope, *(
                        abs(dot(vv, occupancy) - (t-s)/HORIZON)
                        for vv, s, t in zip((vx,vy),state,terminal)
                    ))
                    assert all(abs(t) <= cap+1e-9 for t,cap in zip(terminal,caps))
                    for m,(lo,hi) in enumerate(intervals):
                        covered[m] += lo-1e-9 <= dot(occupancy, METRICS[m]) <= hi+1e-9
                    if a[0] == a[2] and a[1] == a[3]:
                        assert abs(dot(occupancy, METRICS[1]) - a[0]/(a[0]+a[1])) <= 2*caps[0]/HORIZON/(a[0]+a[1])+1e-10
                    if b[0] == b[1] and b[2] == b[3]:
                        assert abs(dot(occupancy, METRICS[2]) - b[0]/(b[0]+b[2])) <= 2*caps[1]/HORIZON/(b[0]+b[2])+1e-10
                    errors["product_11"].append(abs(occupancy[3]-product[3]))
                    errors["entropy_11"].append(abs(occupancy[3]-entropy[3]))
                    errors["unconstrained_interval_11"].append(1.0)
                    for name,(lo,hi) in zip(("balance_11_width","balance_first_width","balance_second_width"),
                                            intervals):
                        errors[name].append(hi-lo)
                    cohort.append(occupancy[3])
                    trials += 1
                ranges.append(max(cohort)-min(cohort))
    assert max_telescope < 1e-10
    return {"trials":trials,"covered":covered,"maximum_balance_residual":max_telescope,
            "matched_kernel_ranges":ranges, "per_run_errors_and_widths":errors}


def result():
    assert sha256(Path(__file__).with_name("SWITCHING_PROTOCOL.md").read_bytes()).hexdigest() == PROTOCOL_HASH
    out = {}
    for index,(label, cohort) in enumerate(models().items()):
        rows = [evaluate(config, index*1000+j) for j, config in enumerate(cohort)]
        metrics = {key: median([value for row in rows for value in row["per_run_errors_and_widths"][key]])
                   for key in rows[0]["per_run_errors_and_widths"]}
        ranges = [value for row in rows for value in row["matched_kernel_ranges"]]
        out[label] = {"fields":len(rows), "trials":sum(r["trials"] for r in rows),
                      "covered_11_first_second":[sum(r["covered"][k] for r in rows) for k in range(3)],
                      "maximum_balance_residual":max(r["maximum_balance_residual"] for r in rows),
                      "median_metrics":metrics,"mean_matched_kernel_range":mean(ranges),
                      "maximum_matched_kernel_range":max(ranges),"field_level":rows}
    primary = out["generated_coupled"]
    return {"protocol_sha256":PROTOCOL_HASH,"cohorts":out,
            "primary_gate_passed": primary["covered_11_first_second"][0] == primary["trials"]
            and primary["median_metrics"]["balance_11_width"] <= .8}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write",type=Path)
    parser.add_argument("--verify",type=Path)
    args = parser.parse_args()
    text = json.dumps(result(),indent=2)+"\n"
    if args.verify:
        if text != args.verify.read_text():
            raise ValueError("Switching results differ")
        print("Verified switching results")
    elif args.write:
        args.write.write_text(text)
    else:
        print(text,end="")


if __name__ == "__main__":
    main()
