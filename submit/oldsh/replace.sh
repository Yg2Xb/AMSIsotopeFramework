#!/usr/bin/env bash
# Patch submit folder (non-recursive) paths and executable name.
# Only files directly under /afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/submit are modified.

set -euo pipefail

SUBMIT_DIR="/afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/submit"

OLD_MC="/afs/cern.ch/work/z/zixuan/logMC"
NEW_MC="/afs/cern.ch/work/z/zixuan/logMC"

OLD_ISS="/afs/cern.ch/work/z/zixuan/logISS"
NEW_ISS="/afs/cern.ch/work/z/zixuan/logISS"

OLD_EXE="SampleProduction"
NEW_EXE="SampleProduction"

echo "[INFO] Target directory: $SUBMIT_DIR"
echo "[INFO] Replace:"
echo "       $OLD_MC  ->  $NEW_MC"
echo "       $OLD_ISS ->  $NEW_ISS"
echo "       $OLD_EXE ->  $NEW_EXE (word boundary)"

shopt -s nullglob

# Process only regular files in the submit directory (no recursion)
for f in "$SUBMIT_DIR"/*; do
  [ -f "$f" ] || continue

  # Backup
  cp -p "$f" "$f.bak"

  # Perform in-place substitutions
  # Note: GNU sed supports \b word boundary. If your sed complains, see fallback below.
  sed -i \
    -e "s#${OLD_MC//#/\\#}#${NEW_MC//#/\\#}#g" \
    -e "s#${OLD_ISS//#/\\#}#${NEW_ISS//#/\\#}#g" \
    -e "s#\\b${OLD_EXE}\\b#${NEW_EXE}#g" \
    "$f"

  echo "[OK] Patched: $(basename "$f")"
done

# Verification summary
echo
echo "[VERIFY] Remaining occurrences (should be empty):"
grep -n --color=auto -e "$OLD_MC" -e "$OLD_ISS" -e '\bmain_exe\b' "$SUBMIT_DIR"/* || echo "[OK] None found."

# Fallback note if sed lacks \b support:
# Uncomment the next block and comment out the \b rule above if needed.
#: <<'FALLBACK'
#for f in "$SUBMIT_DIR"/*; do
#  [ -f "$f" ] || continue
#  # Replace whole word SampleProduction using a broader regex that keeps surrounding chars
#  sed -i -E "s#(^|[^A-Za-z0-9_])${OLD_EXE}([^A-Za-z0-9_]|$)#\1${NEW_EXE}\2#g" "$f"
#done
#FALLBACK

echo "[DONE]"