#!/usr/bin/env python3
"""拉取 MGnify v1 API /runs 端点全部 run（已按 run 去重），用于和本项目语料求重合。

输出 /tmp/mgnify_runs.tsv:  accession  secondary_accession  experiment_type  ena_study  is_private
- accession 多为 SRR/ERR/DRR（run 号，可与我们 run-keyed 数据直接匹配），少量老数据是 SRS/ERS sample 号。
- 并发抓取 + 每页重试；失败页记到 /tmp/mgnify_failed_pages.txt。
代理需在环境里设好 http_proxy/https_proxy。
"""
import json, sys, time, urllib.request, urllib.error
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

BASE = "https://www.ebi.ac.uk/metagenomics/api/v1/runs"
PAGE_SIZE = 250
WORKERS = 12
MAX_RETRY = 6
OUT = "/tmp/mgnify_runs.tsv"
FAILED = "/tmp/mgnify_failed_pages.txt"

def fetch(url):
    req = urllib.request.Request(url, headers={"Accept": "application/json",
                                               "User-Agent": "curl/7.x"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.loads(r.read().decode())

def get_page(page):
    url = f"{BASE}?page_size={PAGE_SIZE}&page={page}"
    last = None
    for attempt in range(MAX_RETRY):
        try:
            d = fetch(url)
            rows = []
            for it in d["data"]:
                a = it["attributes"]
                rows.append((a.get("accession", ""),
                             a.get("secondary-accession", "") or "",
                             a.get("experiment-type", "") or "",
                             a.get("ena-study-accession", "") or "",
                             str(a.get("is-private", ""))))
            return page, rows, None
        except Exception as e:
            last = repr(e)
            time.sleep(2 * (attempt + 1))
    return page, None, last

def main():
    # 拿总数算页数
    meta = fetch(f"{BASE}?page_size=1")["meta"]["pagination"]
    count = meta["count"]
    pages = (count + PAGE_SIZE - 1) // PAGE_SIZE
    sys.stderr.write(f"[start] count={count} pages={pages} workers={WORKERS}\n"); sys.stderr.flush()

    lock = threading.Lock()
    done = [0]
    failed = []
    fout = open(OUT, "w")
    fout.write("accession\tsecondary_accession\texperiment_type\tena_study\tis_private\n")

    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = {ex.submit(get_page, p): p for p in range(1, pages + 1)}
        for fut in as_completed(futs):
            page, rows, err = fut.result()
            with lock:
                done[0] += 1
                if rows is None:
                    failed.append(page)
                    sys.stderr.write(f"[FAIL] page {page}: {err}\n")
                else:
                    for r in rows:
                        fout.write("\t".join(r) + "\n")
                if done[0] % 100 == 0 or done[0] == pages:
                    sys.stderr.write(f"[progress] {done[0]}/{pages} pages, failed={len(failed)}\n")
                    sys.stderr.flush()
                    fout.flush()
    fout.close()
    if failed:
        with open(FAILED, "w") as f:
            f.write("\n".join(map(str, sorted(failed))) + "\n")
        sys.stderr.write(f"[done] {len(failed)} pages failed -> {FAILED}\n")
    else:
        sys.stderr.write("[done] all pages OK\n")
    sys.stderr.flush()

if __name__ == "__main__":
    main()
