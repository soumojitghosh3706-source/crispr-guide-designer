import csv

BBSI = ("GAAGAC", "GTCTTC")


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def make_oligos(guide):
    guide = guide.upper()
    if guide[0] == "G":          # already starts with G: U6 is satisfied
        top = "CACC" + guide
        bottom = "AAAC" + revcomp(guide)
        extra_g = "no"
    else:                        # add a 5' G; matching C goes on the bottom oligo
        top = "CACCG" + guide
        bottom = "AAAC" + revcomp(guide) + "C"
        extra_g = "yes"
    return top, bottom, extra_g


rows = list(csv.DictReader(open("final_guides.csv")))
out = []
for r in rows:
    g = r["guide_seq"].strip().upper()
    top, bottom, extra_g = make_oligos(g)
    notes = []
    if any(site in top or site in bottom for site in BBSI):
        notes.append("BbsI site present - do not use")
    if "TTTT" in g:
        notes.append("TTTT in guide (Pol III terminator)")
    out.append({
        "guide_id": r["guide_id"],
        "gene": r["gene"],
        "guide_seq": g,
        "extra_5prime_G": extra_g,
        "top_oligo_name": r["guide_id"] + "_top",
        "top_oligo_5to3": top,
        "bottom_oligo_name": r["guide_id"] + "_bot",
        "bottom_oligo_5to3": bottom,
        "warnings": "; ".join(notes) or "none",
    })

with open("oligo_order.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
    w.writeheader()
    w.writerows(out)

for o in out:
    print(o["guide_id"], o["top_oligo_5to3"], o["bottom_oligo_5to3"], o["warnings"])
