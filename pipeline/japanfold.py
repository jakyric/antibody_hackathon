"""Minimal JapanFold client (https://japanfold.aiand.com). The API key is read from JAPANFOLD_API_KEY and is
never written to disk or logged."""
from __future__ import annotations

import os
import time
import uuid
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import requests


class JapanFold:
    def __init__(self, base_url: str):
        key = os.environ.get("JAPANFOLD_API_KEY")
        if not key:
            raise SystemExit("Set JAPANFOLD_API_KEY in your shell (get one at https://japanfold.aiand.com/account/token).")
        self.base = base_url.rstrip("/")
        self.s = requests.Session()
        self.s.headers.update({"Authorization": f"Bearer {key}", "X-JapanFold-Client": "her2-pipeline"})

    def submit(self, model: str, name: str, targets: list[dict], params: dict) -> str:
        r = self.s.post(f"{self.base}/v1/predictions", timeout=300, headers={"Idempotency-Key": str(uuid.uuid4())},
                        json={"model": model, "name": name, "params": params, "targets": targets})
        r.raise_for_status()
        return r.json()["id"]

    def wait(self, job: str, poll: int = 55) -> dict:
        while True:
            r = self.s.get(f"{self.base}/v1/jobs/{job}", params={"wait": poll}, timeout=poll + 60).json()
            print(f"    job {job[:8]}: {r.get('status')} {r.get('progress')}")
            if r.get("status") not in ("queued", "running"):
                return r

    def results(self, job: str) -> dict:
        return self.s.get(f"{self.base}/v1/jobs/{job}/results", timeout=120).json()

    def fetch(self, job: str, filename: str, dest: Path, tries: int = 4) -> None:
        if dest.exists() and dest.stat().st_size > 0:
            return
        for t in range(tries):
            try:
                r = self.s.get(f"{self.base}/v1/jobs/{job}/artifacts/{filename}", timeout=300)
                r.raise_for_status()
                dest.write_bytes(r.content)
                return
            except requests.RequestException:
                time.sleep(2 * (t + 1))
        raise RuntimeError(f"could not download {filename}")

    def fetch_many(self, job: str, filenames: list[str], dest_dir: Path, workers: int = 8) -> None:
        dest_dir.mkdir(parents=True, exist_ok=True)
        with ThreadPoolExecutor(workers) as ex:
            list(ex.map(lambda f: self.fetch(job, f, dest_dir / f), filenames))
