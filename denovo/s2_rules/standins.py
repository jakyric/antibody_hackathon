"""Stand-in VHHs used to test the rule harness before real S1 designs exist. NOT designs_*.csv rows."""
import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1]))
import immunum
import dn_common as dn

BCII10 = "QVQLVESGGGLVQPGGSLRLSCAASGGSEYSYSTFSLGWFRQAPGQGLEAVAAIASMGGLTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAAVRGYFMRLPSSHNFRYWGQGTLVTVSS"
TRAS_VH = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS")
A = immunum.Annotator(chains=["h", "k", "l"], scheme="imgt", min_confidence=0.1)


def graft(seq, cdr3=None, cdr1=None, cdr2=None):
    s = A.segment(seq)
    return "".join([s.fr1, cdr1 or s.cdr1, s.fr2, cdr2 or s.cdr2, s.fr3, cdr3 or s.cdr3, s.fr4])


def rows():
    out = []
    def add(i, seq, note):
        s = A.segment(seq)
        out.append(dict(design_id=i, parent_id="", generator="standin", framework="bcii10" if "bcii" in i else "x",
                        seed=0, hotspot_set="", vhh_seq=seq, cdr1=s.cdr1, cdr2=s.cdr2, cdr3=s.cdr3, cdr_scheme="imgt",
                        design_cif="", job_id="", params_json=note, created_utc=dn.utcnow()))
    add("si-bcii10-newcdr3a", graft(BCII10, cdr3="AADLRWGSTYYSDYEYDY"), "BCII10-FGLA framework + invented CDR3")
    add("si-bcii10-newcdr3b", graft(BCII10, cdr3="AAGRYSSGWYPEYDY", cdr1="GFTFSSYAMS"[:10]), "BCII10 FR + invented CDR1/3")
    add("si-bcii10-wt", BCII10, "NbBCII10-FGLA wild type (3EAK), a known antibody")
    add("si-trastuzumab-vh", TRAS_VH, "trastuzumab VH submitted as a vhh")
    add("si-extra-cys", graft(BCII10, cdr3="AADLRWGCTYYSDYEYDY"), "unpaired Cys in CDR3")
    return out

if __name__ == "__main__":
    r = rows()
    dn.write_csv(dn.DENOVO / "s2_rules" / "standins.csv", r, dn.DESIGN_COLS)
    for x in r: print(x["design_id"], len(x["vhh_seq"]), x["cdr1"], x["cdr2"], x["cdr3"])
