"""A thin, idempotent wrapper around forge.py.

Three rules from the platform skill are enforced here rather than left to the
caller, because every one of them costs real time at a hackathon:

1. Never submit the same job twice. Every call carries a label; a label that
   already has a receipt with a run id is *followed*, not resubmitted.
   forge.py does this itself for `submit`, but not for `run`, so we do it for
   both and get one consistent rule.
2. Never guess an input format. Models with a tested recipe are called
   directly; anything else goes through `contract()`, which reads
   input_artifact_contract out of the schema and builds the entry from it.
3. Say the wait before it starts. Every submission prints its estimate.

The run id is read back from forge.py's own receipt file rather than scraped
from stdout, so a change in forge.py's wording does not break us.
"""
import json
import os
import re
import shlex
import subprocess
import sys

from . import config

OP_RE = re.compile(r"\bop-[0-9a-zA-Z]+\b")


class ForgeError(RuntimeError):
    pass


class Forge:
    def __init__(self, forge_py, creds=None, work_dir=None, dry_run=False,
                 budget=1800, verbose=True):
        self.forge_py = os.path.abspath(os.path.expanduser(forge_py))
        self.creds = os.path.abspath(os.path.expanduser(creds)) if creds else None
        self.work_dir = os.path.abspath(work_dir or config.WORK_DIR)
        self.dry_run = dry_run
        self.budget = budget
        self.verbose = verbose
        if not self.dry_run and not os.path.exists(self.forge_py):
            raise ForgeError("forge.py not found at %s" % self.forge_py)
        os.makedirs(os.path.join(self.work_dir, config.INPUTS_DIR), exist_ok=True)
        os.makedirs(os.path.join(self.work_dir, config.RESULTS_DIR), exist_ok=True)

    # ------------------------------------------------------------- plumbing

    def _cmd(self, args):
        cmd = [sys.executable, self.forge_py]
        if self.creds:
            cmd += ["--creds", self.creds]
        return cmd + args

    def _call(self, args, timeout=None):
        cmd = self._cmd(args)
        if self.dry_run:
            print("  [dry-run] cd %s && %s" % (
                shlex.quote(self.work_dir), " ".join(shlex.quote(c) for c in cmd)))
            return ""
        env = dict(os.environ)
        if self.creds:
            env["NEBIUS_CREDS"] = self.creds
        proc = subprocess.run(cmd, cwd=self.work_dir, env=env, text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              timeout=timeout or (self.budget + 300))
        out = proc.stdout or ""
        if self.verbose and out:
            print("\n".join("    " + line for line in out.strip().splitlines()))
        if proc.returncode != 0:
            raise ForgeError("%s failed:\n%s" % (" ".join(args[:2]), out.strip()[-1200:]))
        return out

    # ------------------------------------------------------------- receipts

    def receipt(self, label):
        """forge.py's own record of a run, or None. This is how we avoid
        resubmitting and how `resume` can find the idempotency key later."""
        path = os.path.join(self.work_dir, config.RUNS_DIR, label + ".json")
        if not os.path.exists(path):
            return None
        try:
            with open(path) as fh:
                return json.load(fh)
        except (OSError, json.JSONDecodeError):
            return None

    def existing_op(self, label):
        rec = self.receipt(label) or {}
        return rec.get("operation_id")

    def _op_from(self, label, out):
        op = self.existing_op(label)
        if op:
            return op
        found = OP_RE.search(out or "")
        return found.group(0) if found else None

    # ------------------------------------------------------------ discovery

    def models(self, all_models=False):
        return self._call(["models"] + (["--all"] if all_models else []), timeout=120)

    def schema(self, model, full=False):
        out = self._call(["schema", model] + (["--full"] if full else []), timeout=180)
        try:
            return json.loads(out[out.index("{"):out.rindex("}") + 1])
        except (ValueError, json.JSONDecodeError):
            return {}

    def contract(self, model, operation=None):
        """Derive the --entry shape for a model that has no tested recipe.

        Returns [{name, semantic_type, media_type}, ...] for the operation, as
        published by the model itself. This is the documented way to call a
        model that RECIPES.md does not cover: read the contract, build exactly
        that, submit once.
        """
        schema = self.schema(model, full=True)
        contract = schema.get("input_artifact_contract") or {}
        ops = contract.get("operations") or contract
        if operation and isinstance(ops, dict) and operation in ops:
            ops = ops[operation]
        entries = []
        if isinstance(ops, dict):
            found = ops.get("entries") or ops.get("manifest_entries") or []
            for e in found if isinstance(found, list) else []:
                entries.append({
                    "name": e.get("name"),
                    "semantic_type": e.get("semantic_type") or e.get("type"),
                    "media_type": e.get("media_type"),
                })
        return entries

    # ------------------------------------------------------------ submission

    def run(self, model, args, label, no_wait=False, est=None):
        """Inline model: one call with JSON arguments."""
        if model in config.AVOID:
            raise ForgeError("%s is on the avoid list, see config.AVOID" % model)
        op = self.existing_op(label)
        if op:
            print("  %s: already ran as %s, following it" % (label, op))
            return op
        self._announce(model, label, est)
        path = os.path.join(self.work_dir, config.INPUTS_DIR, label + ".args.json")
        if not self.dry_run:
            with open(path, "w") as fh:
                json.dump(args, fh, indent=1)
        else:
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path, "w") as fh:
                json.dump(args, fh, indent=1)
        rel = os.path.relpath(path, self.work_dir)
        cmd = ["run", model, "--args-file", rel, "--label", label,
               "--budget", str(self.budget)]
        if no_wait:
            cmd.append("--no-wait")
        out = self._retrying(cmd, label)
        return self._op_from(label, out)

    def submit(self, model, entries, params, label, no_wait=False,
               wrap_request=False, est=None):
        """File-based model: upload entries, then submit.

        `entries` is a list of (name, semantic_type, path) with paths relative
        to the working directory.
        """
        if model in config.AVOID:
            raise ForgeError("%s is on the avoid list, see config.AVOID" % model)
        op = self.existing_op(label)
        if op:
            print("  %s: already ran as %s, following it" % (label, op))
            return op
        self._announce(model, label, est)
        cmd = ["submit", model]
        for name, semantic, path in entries:
            cmd += ["--entry", "%s|%s|%s" % (name, semantic, path)]
        cmd += ["--params", json.dumps(params, sort_keys=True),
                "--label", label, "--budget", str(self.budget)]
        if wrap_request:
            cmd.append("--wrap-request")
        if no_wait:
            cmd.append("--no-wait")
        out = self._retrying(cmd, label)
        return self._op_from(label, out)

    def _announce(self, model, label, est):
        secs = est
        if secs is None:
            for _, (name, _lane, s) in config.MODELS.items():
                if name == model:
                    secs = s
                    break
        if secs:
            how = ("about %d seconds" % secs) if secs < 120 else ("about %d minutes" % round(secs / 60.0))
            print("  %s -> %s, first call usually %s (longer when the queue is busy)"
                  % (label, model, how))
        else:
            print("  %s -> %s" % (label, model))

    def _retrying(self, cmd, label):
        """One retry, and only for the one error that is worth retrying.

        A storage hiccup must be resumed, never resubmitted: resubmitting
        starts a second job and we wait twice for the same answer.
        """
        try:
            return self._call(cmd)
        except ForgeError as exc:
            msg = str(exc)
            op = self.existing_op(label)
            if op and ("private storage is temporarily unavailable" in msg
                       or "Storage hiccup" in msg):
                print("  %s: storage hiccup, resuming %s (not resubmitting)" % (label, op))
                return self._call(["resume", op, "--budget", str(self.budget)])
            raise

    # ------------------------------------------------------------- lifecycle

    def wait(self, op):
        if not op:
            return ""
        return self._call(["wait", op, "--budget", str(self.budget)])

    def status(self, op):
        return self._call(["status", op], timeout=120)

    def fetch(self, op, dest):
        """Download a finished run. Returns the destination directory."""
        if not op:
            return None
        full = os.path.join(self.work_dir, dest)
        if not self.dry_run and os.path.isdir(full) and os.listdir(full):
            print("  %s already downloaded to %s" % (op, dest))
            return full
        self._call(["fetch", op, "--dest", dest])
        return full

    def succeeded(self, label):
        rec = self.receipt(label) or {}
        return rec.get("state") == "succeeded"
