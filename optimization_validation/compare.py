#!/usr/bin/env python3
"""Validate a cesar optimization: same output, less memory/CPU?

Usage:    compare.py [reference_cesar] [cesar] [inputs ...]
Defaults: reference_cesar_macos   ../cesar   ../extra/example*.fa

Works on Linux and macOS (peak RSS + CPU come from os.wait4; no /usr/bin/time).
"""

import os
import sys
import re
import glob
import time
import platform
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(SCRIPT_DIR)

# ru_maxrss is in bytes on macOS, kibibytes on Linux.
RSS_TO_BYTES = 1 if platform.system() == "Darwin" else 1024


def natural_key(path):
    """Sort so example2 comes before example10 (numeric, not lexicographic)."""
    name = os.path.basename(path)
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", name)]


def make_runnable(binary, links):
    """cesar finds its profile tables relative to dirname(argv[0]) (the binary's
    directory, not cwd). If the binary does not sit next to extra/, invoke it via
    a short-lived symlink in the repo root. Returns the path to use as argv[0]."""
    binary = os.path.abspath(binary)
    if os.path.dirname(binary) == ROOT:
        return binary
    link = os.path.join(ROOT, ".cmp_%s.%d" % (os.path.basename(binary), os.getpid()))
    if os.path.lexists(link):
        os.remove(link)
    os.symlink(binary, link)
    links.append(link)
    return link


def exit_code(status):
    if os.WIFEXITED(status):
        return os.WEXITSTATUS(status)
    if os.WIFSIGNALED(status):
        return -os.WTERMSIG(status)
    return None


def run_one(argv0, fa, out_path):
    """Run one cesar; return (exit_code, peak_rss_bytes, wall_s, cpu_s)."""
    with open(out_path, "wb") as out, open(os.devnull, "wb") as err:
        start = time.monotonic()
        proc = subprocess.Popen(
            [argv0, fa], stdout=out, stderr=err,
            stdin=subprocess.DEVNULL, cwd=ROOT,
        )
        _, status, ru = os.wait4(proc.pid, 0)
        wall = time.monotonic() - start
        proc.returncode = exit_code(status)  # keep Popen from re-waiting
    return exit_code(status), ru.ru_maxrss * RSS_TO_BYTES, wall, ru.ru_utime + ru.ru_stime


def files_equal(a, b):
    try:
        if os.path.getsize(a) != os.path.getsize(b):
            return False
    except OSError:
        return False
    with open(a, "rb") as fa, open(b, "rb") as fb:
        return fa.read() == fb.read()


def mb(n):
    return "%.1f" % (n / 1048576.0)


def pct(ref, new):
    return "%+.1f%%" % ((new - ref) * 100.0 / ref) if ref > 0 else "n/a"


def main(argv):
    ref = argv[0] if len(argv) > 0 else os.path.join(SCRIPT_DIR, "reference_cesar_macos")
    # ref = argv[0] if len(argv) > 0 else os.path.join(SCRIPT_DIR, "reference_cesar_macos")
    new = argv[1] if len(argv) > 1 else os.path.join(ROOT, "cesar")
    inputs = argv[2:] if len(argv) > 2 else sorted(
        glob.glob(os.path.join(ROOT, "extra", "example*.fa")), key=natural_key)

    for b in (ref, new):
        if not (os.path.isfile(b) and os.access(b, os.X_OK)):
            sys.exit("error: binary not found or not executable: %s" % b)
    if not inputs:
        sys.exit("error: no input .fa files found")

    if platform.system() != "Darwin" and "macos" in os.path.basename(ref):
        sys.stderr.write(
            "WARNING: '%s' is a macOS build but you are on %s; rebuild the "
            "reference for this platform and pass it as arg 1.\n\n"
            % (os.path.basename(ref), platform.system()))

    print("reference: %s" % ref)
    print("new:       %s" % new)
    print("(MB = peak RSS; s = wall seconds; cpu = user+sys CPU seconds)\n")

    hdr = "%-14s %-7s %-6s %9s %9s %8s %8s %8s %8s %8s %8s"
    print(hdr % ("example", "output", "exit", "refMB", "newMB", "mem_d",
                 "ref_s", "new_s", "time_d", "ref_cpu", "new_cpu"))
    print("-" * 106)

    tty = sys.stdout.isatty()
    links = []
    tmp_ref = os.path.join(ROOT, ".cmp_ref.out.%d" % os.getpid())
    tmp_new = os.path.join(ROOT, ".cmp_new.out.%d" % os.getpid())
    n_total = n_equal = n_differ = 0
    try:
        ref_run = make_runnable(ref, links)
        new_run = make_runnable(new, links)
        for fa in inputs:
            if not os.path.exists(fa):
                continue
            name = os.path.basename(fa)
            n_total += 1
            # On a terminal, show which case is running (a big one can take a
            # while); the row printed below overwrites it. Skipped when piped.
            if tty:
                sys.stdout.write("%-14s running...\r" % name)
                sys.stdout.flush()

            rex, rrss, rwall, rcpu = run_one(ref_run, os.path.abspath(fa), tmp_ref)
            nex, nrss, nwall, ncpu = run_one(new_run, os.path.abspath(fa), tmp_new)

            equal = (rex == nex) and files_equal(tmp_ref, tmp_new)
            verdict = "EQUAL" if equal else "DIFFER"
            n_equal += equal
            n_differ += not equal

            print(hdr % (name, verdict, "%s/%s" % (rex, nex),
                         mb(rrss), mb(nrss), pct(rrss, nrss),
                         "%.2f" % rwall, "%.2f" % nwall, pct(rwall, nwall),
                         "%.2f" % rcpu, "%.2f" % ncpu))
    finally:
        for p in links + [tmp_ref, tmp_new]:
            if os.path.lexists(p):
                os.remove(p)

    print("-" * 106)
    print("total: %d   equal: %d   differ: %d" % (n_total, n_equal, n_differ))
    return 1 if n_differ else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
