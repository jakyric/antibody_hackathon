"""Command line entry point.

    # see it without spending anything
    python3 -m denovo --forge ~/.claude/skills/nebius-scientific-ai/forge.py --dry-run

    # the real campaign
    export NEBIUS_CREDS=/path/to/credentials.private.json
    python3 -m denovo --forge ~/.claude/skills/nebius-scientific-ai/forge.py

    # pick up where it stopped (nothing is resubmitted)
    python3 -m denovo --forge ... --stages score,confirm,report
"""
import argparse
import os
import sys

from . import config, report
from .forge_client import Forge, ForgeError
from .pipeline import Campaign

ALL_STAGES = ["target", "msa", "baseline", "generate", "triage", "build",
              "cofold", "score", "confirm", "report"]
DEFAULT_STAGES = ["target", "msa", "baseline", "generate", "triage",
                  "cofold", "score", "confirm", "report"]

DEFAULT_FORGE = "~/.claude/skills/nebius-scientific-ai/forge.py"


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="denovo",
        description="De novo VHH design against %s on the Nebius Scientific AI "
                    "Platform." % config.TARGET_NAME)
    p.add_argument("--forge", default=DEFAULT_FORGE,
                   help="path to forge.py (default: %(default)s)")
    p.add_argument("--creds", default=os.environ.get("NEBIUS_CREDS"),
                   help="credentials.private.json; defaults to $NEBIUS_CREDS. "
                        "Never committed, never printed.")
    p.add_argument("--work-dir", default=os.getcwd(),
                   help="where nebius-runs/ and results/ live (default: cwd). "
                        "Keep it the same between runs or resume cannot find a key.")
    p.add_argument("--route", default="antifold", choices=["antifold", "rfd4", "mosaic"],
                   help="how candidates are generated (default: %(default)s)")
    p.add_argument("--stages", default=",".join(DEFAULT_STAGES),
                   help="comma-separated subset of: " + ",".join(ALL_STAGES))
    p.add_argument("--n-generate", type=int, default=config.N_GENERATE)
    p.add_argument("--keep", type=int, default=config.N_COFOLD,
                   help="candidates that survive triage and get co-folded")
    p.add_argument("--confirm", type=int, default=config.N_CONFIRM,
                   help="top candidates that get a second opinion from OpenFold3")
    p.add_argument("--budget", type=int, default=1800,
                   help="seconds to wait on one run before returning (default: %(default)s)")
    p.add_argument("--target-chain", help="override the target chain in the reference PDB")
    p.add_argument("--heavy-chain", help="override the control's heavy chain")
    p.add_argument("--light-chain", help="override the control's light chain")
    p.add_argument("--dry-run", action="store_true",
                   help="print every forge.py command instead of running it")
    p.add_argument("--quiet", action="store_true", help="hide forge.py output")
    a = p.parse_args(argv)

    stages = [s.strip() for s in a.stages.split(",") if s.strip()]
    unknown = [s for s in stages if s not in ALL_STAGES]
    if unknown:
        p.error("unknown stage(s): %s" % ", ".join(unknown))

    if not a.dry_run and not a.creds:
        p.error("no credentials: pass --creds PATH or set NEBIUS_CREDS "
                "(use --dry-run to see the plan without them)")

    try:
        forge = Forge(a.forge, creds=a.creds, work_dir=a.work_dir,
                      dry_run=a.dry_run, budget=a.budget, verbose=not a.quiet)
        c = Campaign(forge, os.path.abspath(a.work_dir),
                     target_chain=a.target_chain, heavy=a.heavy_chain,
                     light=a.light_chain)

        if "target" in stages:
            c.stage_target()
        if "msa" in stages:
            c.stage_msa()
        if "baseline" in stages:
            c.stage_baseline()
        if "generate" in stages:
            c.stage_generate(route=a.route, n=a.n_generate)
        if "triage" in stages:
            c.stage_triage(keep=a.keep)
        if "build" in stages:
            c.stage_build()
        if "cofold" in stages:
            c.stage_cofold()
        if "score" in stages:
            c.stage_score()
        if "confirm" in stages:
            c.stage_confirm(k=a.confirm)
        if "report" in stages:
            csv_path, md_path = report.write(c.state, os.path.abspath(a.work_dir))
            print("\nwrote %s\nwrote %s" % (csv_path, md_path))
    except ForgeError as exc:
        print("error: %s" % exc, file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        print("\ninterrupted. Nothing was lost: run ids are in %s, state in %s/%s."
              % (config.RUNS_DIR, config.RESULTS_DIR, "campaign.json"), file=sys.stderr)
        return 130
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
