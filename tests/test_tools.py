"""Tests for the helper scripts in tools/: recording runs and scanning before publication."""

import importlib.util
import json
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
TOOLS = ROOT / "tools"


def load(name: str, relative: str):
    spec = importlib.util.spec_from_file_location(name, TOOLS / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


record_run = load("record_run", "record_run.py")
scan_repo = load("scan_repo", "scan_repo.py")


def git(root: Path, *args: str) -> None:
    subprocess.run(["git", "-c", "user.name=T", "-c", "user.email=t@t.no", *args], cwd=root, check=True, capture_output=True)


# ── record_run ──────────────────────────────────────────────────────────────
def test_run_is_recorded_with_seed_inputs_and_ledger(tmp_path):
    (tmp_path / "data.csv").write_text("a,b\n1,2\n")
    run_dir, code = record_run.record(
        "my name!", ["bash", "-c", 'echo "seed=$RUN_SEED" > "$RUN_DIR/out.txt"; echo hello'],
        tmp_path, seed=7, inputs=["data.csv", "missing.csv"], note="trial",
    )
    manifest = json.loads((run_dir / "run.json").read_text())
    assert code == 0 and manifest["exit_code"] == 0
    assert (run_dir / "out.txt").read_text().strip() == "seed=7"
    assert (run_dir / "stdout.log").read_text().strip() == "hello"
    assert manifest["seed"] == 7 and manifest["note"] == "trial"
    assert len(manifest["inputs"]["data.csv"]) == 64
    assert manifest["inputs"]["missing.csv"] is None
    assert run_dir.name.endswith("-my-name")
    rows = (tmp_path / "results" / "ledger.csv").read_text().splitlines()
    assert len(rows) == 2 and rows[1].split(",")[1] == "my name!"


def test_exit_code_is_passed_on_and_ledger_appends(tmp_path):
    _, first = record_run.record("fail", ["bash", "-c", "exit 3"], tmp_path)
    _, second = record_run.record("ok", ["true"], tmp_path)
    assert (first, second) == (3, 0)
    assert len((tmp_path / "results" / "ledger.csv").read_text().splitlines()) == 3


def test_environment_is_not_written_to_the_manifest(tmp_path, monkeypatch):
    monkeypatch.setenv("MY_SECRET", "do-not-write-me")
    run_dir, _ = record_run.record("env", ["true"], tmp_path)
    assert "do-not-write-me" not in (run_dir / "run.json").read_text()


def test_git_state_is_recorded_and_dirty_is_noticed(tmp_path):
    git(tmp_path, "init", "-q")
    (tmp_path / "a.txt").write_text("1")
    git(tmp_path, "add", "a.txt")
    git(tmp_path, "commit", "-qm", "start")
    clean, _ = record_run.record("clean", ["true"], tmp_path)
    (tmp_path / "a.txt").write_text("2")
    dirty, _ = record_run.record("dirty", ["true"], tmp_path)
    assert json.loads((clean / "run.json").read_text())["dirty"] is False
    assert json.loads((dirty / "run.json").read_text())["dirty"] is True


def test_empty_name_is_refused():
    with pytest.raises(ValueError):
        record_run.slugify("!!!")


def test_main_requires_a_command():
    assert record_run.main(["--name", "x"]) == 2
    assert record_run.main(["--name", "x", "--"]) == 2


# ── scan_repo ───────────────────────────────────────────────────────────────
@pytest.fixture
def repo(tmp_path):
    git(tmp_path, "init", "-q")
    (tmp_path / ".gitignore").write_text("*.private.json\nbig.bin\n")
    return tmp_path


def test_clean_repo_has_no_findings(repo):
    (repo / "README.md").write_text("# Hei\n")
    assert scan_repo.scan(repo) == []


def test_ignored_credentials_file_is_not_flagged_but_a_committable_one_is(repo):
    (repo / "credentials.private.json").write_text("{}")
    assert scan_repo.scan(repo) == []
    (repo / "credentials.json").write_text("{}")
    assert any("credentials.json" in f for f in scan_repo.scan(repo))


def test_private_key_and_token_lines_are_found_without_printing_the_value(repo):
    # Assembled here so the test file itself doesn't look like a leak to the scanner.
    fake_key = "-----BEGIN " + "OPENSSH PRIVATE KEY-----"
    fake_field = "api_" + "key = 'abcdefghijklmnop1234'"
    (repo / "notes.txt").write_text(f"{fake_key}\n{fake_field}\n")
    findings = scan_repo.scan(repo)
    assert any("notes.txt:1" in f and "private key" in f for f in findings)
    assert any("notes.txt:2" in f for f in findings)
    assert not any("abcdefghijklmnop1234" in f for f in findings)


def test_env_files_and_keys_are_flagged_by_name(repo):
    (repo / ".env").write_text("A=1")
    (repo / "id_ed25519").write_text("x")
    names = " ".join(scan_repo.scan(repo))
    assert ".env" in names and "id_ed25519" in names


def test_large_files_are_flagged_unless_ignored(repo):
    for name in ("big.bin", "big2.bin"):
        with (repo / name).open("wb") as handle:
            handle.truncate(51 * 1024 * 1024)
    findings = scan_repo.scan(repo)
    assert any("big2.bin" in f for f in findings)
    assert not any(f.startswith("big.bin") for f in findings)
