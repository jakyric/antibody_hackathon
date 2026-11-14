"""Try to pull HER2 VHH sequences printed in the open-access papers listed in research/denovo_vhh_research.md (Q1).
Only sequences literally present in the paper text are kept. Writes cache/lit_candidates.json for manual review."""
import json, re, urllib.request, urllib.parse, sys
PMIDS = {"37632647": "Nikkhoi2023", "30518413": "Hussack2018", "22440788": "RahimiJamnani2012", "37954614": "Liu2023",
         "33738081": "Kim2021"}
SEQ = re.compile(r"[EQD][VIL]QL[A-Z]{95,135}?[VT]SS")
def get(u):
    return urllib.request.urlopen(urllib.request.Request(u, headers={"User-Agent": "s2-rules/0.1"}), timeout=60).read().decode(errors="replace")
out = {}
for pmid, tag in PMIDS.items():
    try:
        j = json.loads(get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:{pmid}%20AND%20SRC:MED&format=json&resultType=core"))
        rec = j["resultList"]["result"][0]
        pmcid = rec.get("pmcid")
        oa = rec.get("isOpenAccess")
        print(tag, pmid, pmcid, "OA" if oa == "Y" else "not-OA", rec.get("title", "")[:70])
        if not pmcid:
            continue
        xml = get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML")
        txt = re.sub(r"<[^>]+>", " ", xml)
        txt = re.sub(r"\s+", "", txt)  # sequences may be split by whitespace/line breaks in tables
        seqs = sorted(set(m.group(0) for m in SEQ.finditer(txt)))
        print("   candidate VHH-like sequences in text:", len(seqs))
        out[tag] = dict(pmid=pmid, pmcid=pmcid, seqs=seqs)
    except Exception as e:
        print(tag, pmid, "failed:", str(e)[:100])
json.dump(out, open("cache/lit_candidates.json", "w"), indent=1)
