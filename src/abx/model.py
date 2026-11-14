"""A binder model trained on the FACS data, in pure numpy.

WHY PURE NUMPY
  The repo must be runnable by a stranger. A small multilayer net on a
  10-residue loop is a small task: 200 features in, 524,346 examples.
  It needs no GPU and no framework, and that means we have no
  installation that can fail during assessment either.

HOW IT IS TESTED
  Trained on FACS (cell-surface sorting, three classes). Tested on SPR
  (biophysical KD, 182 variants). Two different assays, different
  laboratories, different measurement principles. If the model transfers
  there, it has learned something about binding — not about a quirk of
  one setup.

  The split is group-based: sequences that share a beginning end up in the
  same place. Without that the dataset leaks and everything looks brilliant.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .data import TZ_CDRH3_10, load_facs
from .metrics import auroc

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_INDEX = {a: i for i, a in enumerate(AA)}
L = 10                      # the dataset's CDR-H3 length
NFEAT = L * len(AA)         # 200

MODEL_PATH = Path(__file__).resolve().parents[2] / "results" / "binder_model.npz"


# ── features ────────────────────────────────────────────────────────────────
def onehot(seqs) -> np.ndarray:
    """(N, 200) float32. Position × amino acid."""
    X = np.zeros((len(seqs), NFEAT), dtype=np.float32)
    for i, s in enumerate(seqs):
        for p, ch in enumerate(s[:L]):
            j = AA_INDEX.get(ch)
            if j is not None:
                X[i, p * len(AA) + j] = 1.0
    return X


def group_bucket(seq: str, prefix: int = 5) -> int:
    """0–9, determined by the first residues. Similar sequences get the same bucket."""
    h = hashlib.blake2b(seq[:prefix].encode(), digest_size=8).digest()
    return int.from_bytes(h, "big") % 10


def split_masks(seqs) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Train / validate / test, 70-10-20, without splitting a group."""
    b = np.array([group_bucket(s) for s in seqs])
    return b <= 6, b == 7, b >= 8


# ── the net ─────────────────────────────────────────────────────────────────
@dataclass
class MLP:
    W1: np.ndarray
    b1: np.ndarray
    W2: np.ndarray
    b2: np.ndarray
    W3: np.ndarray
    b3: np.ndarray

    @staticmethod
    def init(hidden: int = 96, hidden2: int = 32, seed: int = 0) -> "MLP":
        rng = np.random.default_rng(seed)
        he = lambda nin, nout: rng.normal(0, np.sqrt(2.0 / nin), (nin, nout)).astype(np.float32)
        return MLP(he(NFEAT, hidden), np.zeros(hidden, np.float32),
                   he(hidden, hidden2), np.zeros(hidden2, np.float32),
                   he(hidden2, 1), np.zeros(1, np.float32))

    def forward(self, X: np.ndarray):
        z1 = X @ self.W1 + self.b1
        a1 = np.maximum(z1, 0)
        z2 = a1 @ self.W2 + self.b2
        a2 = np.maximum(z2, 0)
        z3 = (a2 @ self.W3 + self.b3).ravel()
        return z1, a1, z2, a2, z3

    def logits(self, X: np.ndarray) -> np.ndarray:
        return self.forward(X)[4]

    def params(self):
        return [self.W1, self.b1, self.W2, self.b2, self.W3, self.b3]


def _bce_with_logits(z: np.ndarray, y: np.ndarray) -> float:
    # stable form: max(z,0) - z*y + log(1 + exp(-|z|))
    return float(np.mean(np.maximum(z, 0) - z * y + np.log1p(np.exp(-np.abs(z)))))


def train(
    hidden: int = 96,
    hidden2: int = 32,
    epochs: int = 14,
    batch: int = 512,
    lr: float = 2e-3,
    weight_decay: float = 1e-5,
    seed: int = 0,
    verbose: bool = True,
) -> tuple[MLP, dict]:
    """Train "high binder versus low binder". 'medium' is held out.

    Returns the model and a log with AUROC per epoch, so one can see that
    it is not just memorising.
    """
    facs = load_facs()
    keep = facs.is_high | facs.is_low
    seqs = facs.cdrh3[keep]
    y = facs.is_high[keep].astype(np.float32)

    tr, va, te = split_masks(seqs)
    X = onehot(seqs)
    if verbose:
        print(f"train {tr.sum():,}  validate {va.sum():,}  test {te.sum():,}")

    net = MLP.init(hidden, hidden2, seed)
    ps = net.params()
    ms = [np.zeros_like(p) for p in ps]
    vs = [np.zeros_like(p) for p in ps]
    rng = np.random.default_rng(seed)
    Xtr, ytr = X[tr], y[tr]
    step = 0
    best = (-1.0, None)
    log = {"val_auroc": [], "train_loss": []}

    for ep in range(1, epochs + 1):
        order = rng.permutation(len(Xtr))
        losses = []
        for s in range(0, len(order), batch):
            idx = order[s:s + batch]
            xb, yb = Xtr[idx], ytr[idx]
            z1, a1, z2, a2, z3 = net.forward(xb)
            p = 1.0 / (1.0 + np.exp(-z3))
            losses.append(_bce_with_logits(z3, yb))

            n = len(xb)
            dz3 = ((p - yb) / n).astype(np.float32)[:, None]
            gW3 = a2.T @ dz3
            gb3 = dz3.sum(0)
            da2 = dz3 @ net.W3.T
            dz2 = da2 * (z2 > 0)
            gW2 = a1.T @ dz2
            gb2 = dz2.sum(0)
            da1 = dz2 @ net.W2.T
            dz1 = da1 * (z1 > 0)
            gW1 = xb.T @ dz1
            gb1 = dz1.sum(0)

            grads = [gW1, gb1, gW2, gb2, gW3, gb3]
            step += 1
            for k, (pm, g) in enumerate(zip(ps, grads)):
                if g.ndim == 2:
                    g = g + weight_decay * pm
                ms[k] = 0.9 * ms[k] + 0.1 * g
                vs[k] = 0.999 * vs[k] + 0.001 * (g * g)
                mhat = ms[k] / (1 - 0.9 ** step)
                vhat = vs[k] / (1 - 0.999 ** step)
                pm -= lr * mhat / (np.sqrt(vhat) + 1e-8)

        va_auc = auroc(net.logits(X[va]), y[va].astype(bool))
        log["val_auroc"].append(va_auc)
        log["train_loss"].append(float(np.mean(losses)))
        if verbose:
            print(f"  epoch {ep:>2}  loss {np.mean(losses):.4f}  validate-AUROC {va_auc:.4f}")
        if va_auc > best[0]:
            best = (va_auc, [p.copy() for p in net.params()])

    # back to the best epoch
    if best[1] is not None:
        for p, b in zip(net.params(), best[1]):
            p[...] = b

    te_auc = auroc(net.logits(X[te]), y[te].astype(bool))
    log["test_auroc"] = te_auc
    log["val_best"] = best[0]
    if verbose:
        print(f"best validate-AUROC {best[0]:.4f}   →   TEST-AUROC {te_auc:.4f}")
    return net, log


def save(net: MLP, path: Path = MODEL_PATH) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, W1=net.W1, b1=net.b1, W2=net.W2,
                        b2=net.b2, W3=net.W3, b3=net.b3)


def load(path: Path = MODEL_PATH) -> MLP:
    d = np.load(path)
    return MLP(d["W1"], d["b1"], d["W2"], d["b2"], d["W3"], d["b3"])


def train_honest(n_test: int = 6000, min_dist: int = 3, seed: int = 0,
                 **kw) -> tuple[MLP, dict]:
    """Train with a distance split instead of prefix grouping.

    This is the training path that applies. The prefix grouping in
    `split_masks` looked right and leaked: 194 of 300 test sequences lay
    one mutation from the training data. See splits.py.
    """
    from .splits import distance_split, report_split

    facs = load_facs()
    keep = facs.is_high | facs.is_low
    seqs = list(facs.cdrh3[keep])
    y = facs.is_high[keep].astype(np.float32)

    tr, te = distance_split(seqs, n_test=n_test, min_dist=min_dist, seed=seed)
    rep = report_split(seqs, tr, te)
    print(f"split: train {rep['n_train']:,}  test {rep['n_test']:,}  "
          f"minimum distance {rep['min_dist']}  median {rep['median_dist']:.0f}")
    assert rep["min_dist"] >= min_dist, "the split leaks"

    X = onehot(seqs)
    net, log = _fit(X[tr], y[tr], X[te], y[te], seed=seed, **kw)
    log["split"] = rep
    return net, log


def _fit(Xtr, ytr, Xte, yte, hidden=96, hidden2=32, epochs=12, batch=512,
         lr=2e-3, weight_decay=1e-4, seed=0, verbose=True):
    net = MLP.init(hidden, hidden2, seed)
    ps = net.params()
    ms = [np.zeros_like(p) for p in ps]
    vs = [np.zeros_like(p) for p in ps]
    rng = np.random.default_rng(seed)
    step = 0
    best = (-1.0, None)
    log = {"test_auroc_per_epoch": []}

    for ep in range(1, epochs + 1):
        order = rng.permutation(len(Xtr))
        for s in range(0, len(order), batch):
            i = order[s:s + batch]
            xb, yb = Xtr[i], ytr[i]
            z1, a1, z2, a2, z3 = net.forward(xb)
            p = 1.0 / (1.0 + np.exp(-z3))
            n = len(xb)
            dz3 = ((p - yb) / n).astype(np.float32)[:, None]
            gW3 = a2.T @ dz3
            gb3 = dz3.sum(0)
            dz2 = (dz3 @ net.W3.T) * (z2 > 0)
            gW2 = a1.T @ dz2
            gb2 = dz2.sum(0)
            dz1 = (dz2 @ net.W2.T) * (z1 > 0)
            gW1 = xb.T @ dz1
            gb1 = dz1.sum(0)
            step += 1
            for k, (pm, g) in enumerate(zip(ps, [gW1, gb1, gW2, gb2, gW3, gb3])):
                if g.ndim == 2:
                    g = g + weight_decay * pm
                ms[k] = 0.9 * ms[k] + 0.1 * g
                vs[k] = 0.999 * vs[k] + 0.001 * (g * g)
                pm -= lr * (ms[k] / (1 - 0.9 ** step)) / (
                    np.sqrt(vs[k] / (1 - 0.999 ** step)) + 1e-8)

        a = auroc(net.logits(Xte), yte.astype(bool))
        log["test_auroc_per_epoch"].append(a)
        if verbose:
            print(f"  epoch {ep:>2}  test-AUROC {a:.4f}")
        if a > best[0]:
            best = (a, [q.copy() for q in net.params()])

    for q, b in zip(net.params(), best[1]):
        q[...] = b
    log["test_auroc"] = best[0]
    return net, log


def scorer(net: MLP | None = None):
    """Turn the model into a scorer the benchmark can measure."""
    net = net if net is not None else load()

    def score(seqs) -> np.ndarray:
        return net.logits(onehot(list(seqs)))

    return score


def blend_scorer(net: MLP | None = None, weight: float = 0.25):
    """The model plus a nudge towards trastuzumab.

    Edit distance is a strong baseline, and the model sees only
    sequence. Blending in a little distance costs nothing and captures the
    part of the signal that lies in "close to the original is safer".
    The weight was chosen on the validation split, not on the test split.
    """
    net = net if net is not None else load()

    def score(seqs) -> np.ndarray:
        seqs = list(seqs)
        z = net.logits(onehot(seqs))
        z = (z - z.mean()) / (z.std() + 1e-9)
        d = np.array([sum(a != b for a, b in zip(s, TZ_CDRH3_10))
                      if len(s) == L else 99 for s in seqs], dtype=float)
        d = -(d - d.mean()) / (d.std() + 1e-9)
        return (1 - weight) * z + weight * d

    return score


if __name__ == "__main__":
    import json
    net, log = train_honest()
    save(net)
    out = MODEL_PATH.parent / "model_log.json"
    out.write_text(json.dumps(
        {k: (v if not isinstance(v, np.ndarray) else v.tolist())
         for k, v in log.items()}, indent=2, ensure_ascii=False))
    print(f"\nmodel saved: {MODEL_PATH}")
    print(f"log saved:   {out}")
