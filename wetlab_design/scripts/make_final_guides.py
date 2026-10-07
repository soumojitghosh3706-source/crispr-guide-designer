import csv

FILES = {
    "PIP4K2C": "phase3_comparison.csv",
    "MUTYH": "mutyh_phase3_comparison.csv",
    "KRAS": "kras_phase3_comparison.csv",
    "G6PD": "g6pd_phase3_comparison.csv",
}
COLS = ["gene", "guide_id", "guide_seq", "exon", "our_verdict",
        "mit_specificity", "cfd_specificity", "flag", "variant"]

rows = []
for gene, path in FILES.items():
    n = 0
    for r in csv.DictReader(open(path)):
        if r["comparison"].strip() == "Agreement":
            n += 1
            rows.append({
                "gene": gene,
                "guide_id": f"{gene}_g{n}",
                "guide_seq": r["guide"].strip().upper(),
                "exon": r["exon"],
                "our_verdict": r["our_verdict"],
                "mit_specificity": r["mit_specificity"],
                "cfd_specificity": r["cfd_specificity"],
                "flag": r.get("flag", ""),
                "variant": r.get("variant", ""),
            })

with open("final_guides.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=COLS)
    w.writeheader()
    w.writerows(rows)
print(len(rows), "final guides written")
