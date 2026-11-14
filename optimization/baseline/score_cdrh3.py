"""Gate: score trastuzumab CDRH3 candidates for HER2 binding likelihood and the novelty rule.

  optimization/baseline/.venv/bin/python optimization/baseline/score_cdrh3.py --seqs WGGDGFYAMD WGSQGFYAYD
  optimization/baseline/.venv/bin/python optimization/baseline/score_cdrh3.py --csv designs.csv --col cdrh3 --out scored.csv
  optimization/baseline/.venv/bin/python optimization/baseline/score_cdrh3.py --fasta cdrh3.fa

Output columns
  p_hgb         P(high HER2 binder) from our gradient-boosted model (trained on all FACS variants)
  p_authors_cnn P(high) from the Mehta/Chinery CNN (their weights, numpy re-implementation)
  p_mean        average of the two
  n_changed     CDRH3 positions that differ from trastuzumab (WGGDGFYAMD)
  identity_wt   1 - n_changed / 10
  near_provided True if a FACS or SPR sequence is within 1 mutation (the rule needs >= 2 away)
  challenge_ok  identity_wt <= 0.70 and not near_provided  (Apheris docs' "challenge filter";
                abguard remains the official checker)
  challenge_ok_strict  challenge_ok and identity_wt < 0.70 (task PDF wording: "under 70%", so >= 4 changes)
  keep          challenge_ok_strict and p_mean >= --threshold

Only length-10 CDRH3 strings are supported. The score says "likely a high binder in the yeast
display screen"; it is NOT an affinity predictor (see baseline/README.md).
"""
import argparse
import sys

import joblib
import numpy as np
import pandas as pd

import common as C


def read_inputs(a):
    if a.seqs:
        return pd.DataFrame({"cdrh3": a.seqs})
    if a.csv:
        df = pd.read_csv(a.csv)
        return df.rename(columns={a.col: "cdrh3"})
    if a.fasta:
        recs, name, buf = [], None, []
        for line in open(a.fasta):
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                name, buf = line[1:], []
            elif line:
                buf.append(line)
        if name is not None:
            recs.append((name, "".join(buf)))
        return pd.DataFrame(recs, columns=["id", "cdrh3"])
    sys.exit("give --seqs, --csv or --fasta")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seqs", nargs="+")
    ap.add_argument("--csv")
    ap.add_argument("--col", default="cdrh3")
    ap.add_argument("--fasta")
    ap.add_argument("--out")
    ap.add_argument("--threshold", type=float, default=0.5)
    ap.add_argument("--skip-novelty", action="store_true", help="skip loading FACS/SPR for the novelty check")
    a = ap.parse_args()

    df = read_inputs(a)
    df["cdrh3"] = df["cdrh3"].astype(str).str.strip().str.upper()
    bad = df.cdrh3.str.len() != C.L
    if bad.any():
        print(f"warning: dropping {int(bad.sum())} sequence(s) that are not length {C.L}", file=sys.stderr)
        df = df[~bad].reset_index(drop=True)
    idx = C.encode_idx(df.cdrh3)

    # our own model, written by train_eval.py (never load pickles from other people)
    bundle = joblib.load(C.OUT / "hgb_full.joblib")
    df["p_hgb"] = bundle["model"].predict_proba(idx.astype(np.float32))[:, 1]
    df["p_authors_cnn"] = C.AuthorsCNN(C.MEHTA / "saved_weights.h5").predict_proba(idx)
    df["p_mean"] = (df.p_hgb + df.p_authors_cnn) / 2
    df["n_changed"] = C.hamming_to_wt(idx)
    df["identity_wt"] = 1 - df.n_changed / C.L
    if a.skip_novelty:
        df["near_provided"] = np.nan
        df["challenge_ok"] = df.identity_wt <= 0.70
    else:
        df["near_provided"] = C.ProvidedData(C.provided_sequences()).has_neighbour_within_1(idx)
        df["challenge_ok"] = (df.identity_wt <= 0.70) & ~df.near_provided
    # the task PDF says "under 70%" (strict, so at least 4 of 10 changed); the Apheris docs say "at most 70%"
    df["challenge_ok_strict"] = df.challenge_ok & (df.identity_wt < 0.70)
    df["keep"] = df.challenge_ok_strict & (df.p_mean >= a.threshold)

    cols = [c for c in df.columns if c not in ("p_hgb", "p_authors_cnn", "p_mean")]
    df = df[cols[:1] + ["p_hgb", "p_authors_cnn", "p_mean"] + cols[1:]]
    text = df.round(4).to_csv(index=False)
    if a.out:
        open(a.out, "w").write(text)
        print(f"wrote {len(df)} rows to {a.out}; {int(df.keep.sum())} pass the gate", file=sys.stderr)
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
