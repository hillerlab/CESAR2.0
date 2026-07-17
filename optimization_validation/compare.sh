#!/usr/bin/env bash
#
# Validate a cesar optimization: same output, less memory/CPU?
# Usage:    compare.sh <reference_cesar> <cesar> [inputs...]
# Defaults: reference_cesar_macos   ./cesar   extra/example*.fa

set -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$ROOT" || exit 1

REF="${1:-$SCRIPT_DIR/reference_cesar_macos}"
NEW="${2:-$ROOT/cesar}"
shift 2 2>/dev/null || true
if [ "$#" -gt 0 ]; then INPUTS=("$@"); else INPUTS=("$ROOT"/extra/example*.fa); fi

for bin in "$REF" "$NEW"; do
  if [ ! -x "$bin" ]; then echo "error: binary not found or not executable: $bin" >&2; exit 1; fi
done

# Warn if the reference is a macOS build but we are not on macOS (won't exec).
if [ "$(uname -s)" != "Darwin" ] && [[ "$(basename "$REF")" == *macos* ]]; then
  echo "WARNING: '$(basename "$REF")' is a macOS build but you are on $(uname -s); rebuild the reference for this platform and pass it as arg 1." >&2
  echo >&2
fi

# /usr/bin/time: macOS uses -l and reports RSS in bytes; GNU (Linux) uses -v
# and reports RSS in kbytes. Timing lines also differ; both are parsed below.
TIME_BIN=/usr/bin/time
case "$(uname -s)" in
  Darwin) TIME_FLAG="-l"; RSS_UNIT=1 ;;
  *)      TIME_FLAG="-v"; RSS_UNIT=1024 ;;
esac

tmp="$(mktemp -d)"
LINKS=()
cleanup() { rm -rf "$tmp"; for l in "${LINKS[@]:-}"; do [ -n "$l" ] && rm -f "$l"; done; }
trap cleanup EXIT

# cesar resolves its tables relative to dirname(argv[0]) (the binary's dir, not
# cwd), so a binary must sit next to extra/. If it lives elsewhere, invoke it via
# a short-lived symlink in the repo root. Sets REPLY / appends to LINKS; must run
# in the main shell (not a $(...) subshell) or cleanup would miss the symlink.
make_runnable() {
  local target="$1" link
  if [ "$(dirname "$target")" = "$ROOT" ]; then REPLY="$target"; return; fi
  link="$ROOT/.cmp_$(basename "$target").$$"
  ln -sf "$target" "$link"
  LINKS+=("$link")
  REPLY="$link"
}
make_runnable "$REF"; REF_RUN="$REPLY"
make_runnable "$NEW"; NEW_RUN="$REPLY"

peak_bytes() {  # peak RSS -> bytes (macOS: field is bytes; GNU: kbytes)
  awk -v u="$RSS_UNIT" '/[Mm]aximum resident set size/{
      for (i=1;i<=NF;i++) if ($i ~ /^[0-9]+$/) { print $i*u; exit }
  }' "$1"
}
wall_seconds() {  # elapsed wall-clock seconds (macOS "N real" | GNU "Elapsed ... h:mm:ss")
  awk '
    { for (i=1;i<=NF;i++) if ($i=="real") wm=$(i-1) }
    /Elapsed \(wall clock\)/ { t=$NF }
    END {
      if (wm!="") { print wm+0 }
      else {
        n=split(t,a,":");
        if (n==3) print a[1]*3600+a[2]*60+a[3];
        else if (n==2) print a[1]*60+a[2];
        else print a[1]+0;
      }
    }' "$1"
}
cpu_seconds() {  # user+sys CPU seconds (macOS "N user N sys" | GNU "User/System time")
  awk '
    { for (i=1;i<=NF;i++) { if ($i=="user") um=$(i-1); if ($i=="sys") sm=$(i-1) } }
    /User time \(seconds\)/   { ul=$NF }
    /System time \(seconds\)/ { sl=$NF }
    END { if (um!="" || sm!="") print um+0 + sm+0; else print ul+0 + sl+0 }' "$1"
}
mb()  { awk -v b="${1:-0}" 'BEGIN{ printf "%.1f", b/1048576 }'; }
pct() { awk -v r="${1:-0}" -v n="${2:-0}" 'BEGIN{ if(r>0) printf "%+.1f%%", (n-r)*100.0/r; else printf "n/a" }'; }

echo "reference: $REF"
echo "new:       $NEW"
echo "(MB = peak RSS; s = wall seconds; cpu = user+sys CPU seconds)"
echo
HDR="%-13s %-7s %-6s %8s %8s %8s %8s %8s %8s %8s %8s\n"
printf "$HDR" "example" "output" "exit" \
       "refMB" "newMB" "mem_d" "ref_s" "new_s" "time_d" "ref_cpu" "new_cpu"
printf '%.0s-' {1..104}; echo

n_total=0; n_equal=0; n_differ=0
for fa in "${INPUTS[@]}"; do
  [ -e "$fa" ] || continue
  name="$(basename "$fa")"
  n_total=$((n_total+1))

  $TIME_BIN $TIME_FLAG "$REF_RUN" "$fa" >"$tmp/ref.out" 2>"$tmp/ref.time"; rex=$?
  $TIME_BIN $TIME_FLAG "$NEW_RUN" "$fa" >"$tmp/new.out" 2>"$tmp/new.time"; nex=$?

  if [ "$rex" = "$nex" ] && diff -q "$tmp/ref.out" "$tmp/new.out" >/dev/null 2>&1; then
    verdict="EQUAL"; n_equal=$((n_equal+1))
  else
    verdict="DIFFER"; n_differ=$((n_differ+1))
  fi

  rp="$(peak_bytes  "$tmp/ref.time")"; rp="${rp:-0}"
  np="$(peak_bytes  "$tmp/new.time")"; np="${np:-0}"
  rw="$(wall_seconds "$tmp/ref.time")"; nw="$(wall_seconds "$tmp/new.time")"
  rc="$(cpu_seconds  "$tmp/ref.time")"; nc="$(cpu_seconds  "$tmp/new.time")"

  printf "$HDR" "$name" "$verdict" "$rex/$nex" \
         "$(mb "$rp")" "$(mb "$np")" "$(pct "$rp" "$np")" \
         "$rw" "$nw" "$(pct "$rw" "$nw")" "$rc" "$nc"
done

printf '%.0s-' {1..104}; echo
echo "total: $n_total   equal: $n_equal   differ: $n_differ"
[ "$n_differ" -eq 0 ]
