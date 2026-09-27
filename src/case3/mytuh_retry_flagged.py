"""
mutyh_retry_flagged.py  --  helper for Case study 3, Step 4

Removes specific guides from mutyh_blast_progress.csv (and their rows from
mutyh_blast_hits_genomic.csv, if any) so that re-running
mutyh_blast_offtargets.py retries them from scratch, instead of skipping
them as already "done". Use this to re-verify any guide whose earlier run
might have hit a silent fetch error.

Edit RETRY_GUIDES below, then run this once, then run
mutyh_blast_offtargets.py again as normal.
"""

import os
import pandas as pd

RETRY_GUIDES = {
    "ATCAGGTGGTAGAGGAGCTA",
    "CACAGGAGGTGAATCAACTC",
    "TCACCACACTCCTCCACGTC",
}

PROGRESS_FILE = "mutyh_blast_progress.csv"
HITS_FILE = "mutyh_blast_hits_genomic.csv"

if __name__ == "__main__":
    print(f"Working folder: {os.getcwd()}")
    print(f"Files here: {sorted(f for f in os.listdir('.') if f.startswith('mutyh'))}\n")

    if not os.path.exists(PROGRESS_FILE):
        raise SystemExit(
            f"{PROGRESS_FILE} not found in this folder. This script must run "
            f"from the EXACT same folder as mutyh_blast_offtargets.py -- move "
            f"or copy this file there and run it again."
        )

    p_before = pd.read_csv(PROGRESS_FILE)
    already_done = set(p_before.loc[p_before["status"] == "done", "guide"])
    matched = RETRY_GUIDES & already_done
    print(f"Guides currently marked 'done' that match RETRY_GUIDES: "
         f"{len(matched)} of {len(RETRY_GUIDES)}")
    if len(matched) < len(RETRY_GUIDES):
        print("  Missing:", RETRY_GUIDES - already_done,
             "-- check for typos or already-retried guides.")

    p = p_before[~p_before["guide"].isin(RETRY_GUIDES)]
    p.to_csv(PROGRESS_FILE, index=False)
    print(f"\nProgress rows: {len(p_before)} -> {len(p)} "
         f"(removed {len(p_before) - len(p)})")

    if os.path.exists(HITS_FILE):
        h_before = pd.read_csv(HITS_FILE)
        h = h_before[~h_before["guide"].isin(RETRY_GUIDES)]
        h.to_csv(HITS_FILE, index=False)
        print(f"Hit rows: {len(h_before)} -> {len(h)} "
             f"(removed {len(h_before) - len(h)})")
    else:
        print(f"{HITS_FILE} not found here (fine if no hits were ever saved).")

    # Prove it actually took, by re-reading the file fresh from disk
    p_check = pd.read_csv(PROGRESS_FILE)
    still_there = RETRY_GUIDES & set(p_check["guide"])
    if still_there:
        print(f"\nWARNING: these guides are STILL in {PROGRESS_FILE} after "
             f"saving: {still_there}. Something is wrong -- do not run "
             f"mutyh_blast_offtargets.py yet.")
    else:
        print(f"\nConfirmed: none of the {len(RETRY_GUIDES)} retry guides "
             f"remain in {PROGRESS_FILE}.")
        print("Now run mutyh_blast_offtargets.py again from THIS SAME FOLDER "
             "-- it will retry only these guides and skip the rest.")
             