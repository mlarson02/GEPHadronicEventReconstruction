#!/usr/bin/env bash
# Sweep every memory-print .dat under MemPrints_v3 and write a GEP-format copy of each
# into a parallel MemPrints_v3_GEPFormat tree, preserving the subdirectory layout.
#
# The conversion itself is convert_dat.py: it pulls the trailing 0x word off each object
# line, zero-pads it to the target width, and writes one space-separated line per event.
# This script only decides, per source directory, WHICH width and how many words a line
# may carry -- and those two numbers are properties of the consuming firmware block, not
# of the .dat, so they cannot be inferred from the file and live in the table below.
#
# Both the converted .txt and a copy of the original .dat land in the destination, so the
# new tree is self-contained and a consumer never has to reach back into MemPrints_v3.
#
# Usage:
#   ./convert_all_dat.sh              # convert everything
#   ./convert_all_dat.sh --dry-run    # print what would be done, touch nothing
set -uo pipefail

SRC="/eos/home-m/mlarson/TransferMemPrintsLUTs/data/MemPrints_v3"
DST="/eos/home-m/mlarson/TransferMemPrintsLUTs/data/MemPrints_v3_GEPFormat"
CONVERTER="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/convert_dat.py"

DRY_RUN=0
[ "${1:-}" = "--dry-run" ] && DRY_RUN=1

if [ ! -f "$CONVERTER" ]; then
    echo "Error: converter not found at $CONVERTER" >&2
    exit 1
fi
if [ ! -d "$SRC" ]; then
    echo "Error: source tree not found at $SRC" >&2
    exit 1
fi

# Word width and per-event word cap, by source directory.
#
# maxWords 0 means "no truncation" -- write every word the event has.
#
# JET1 is the only 128-bit type; everything else is 64. The source words are narrower
# than both (the tower and cone-jet prints carry 32-bit words), so the padding in
# convert_dat.py is what brings them up to the target width.
settings_for() {
    case "$1" in
        # JET1 (WTA cone jets): the one 128-bit type, 10 jets per event.
        GEPConeJetsCellsTowers|GEPConeJetsCellsTowersEtaSK)
            bits=128; maxWords=10 ;;
        # Cell towers, including the sorted variant.
        GEPCellsTowers|GEPCellsTowersEtaSK|GEPCellsTowers_Sort)
            bits=64;  maxWords=256 ;;
        # Jet tagger large-R jets: leading and subleading only.
        largeRJets_*)
            bits=64;  maxWords=2 ;;
        # MET: a single TOB per event.
        METBaselineValidation)
            bits=64;  maxWords=1 ;;
        # FEX small-R jets. 64-bit like everything non-JET1; the per-event cap was never
        # specified, so nothing is truncated. Set a number here once it is known.
        gFex|jFex)
            bits=64;  maxWords=0 ;;
        # Anything added to MemPrints_v3 later lands here. 64-bit, untruncated, and
        # announced, so a new directory is noticed rather than silently guessed at.
        *)
            bits=64;  maxWords=0
            echo "  [warn] no rule for directory '$1' -- defaulting to ${bits}-bit, untruncated" >&2 ;;
    esac
}

nConverted=0
nFailed=0

# -print0 / read -d '' so paths with spaces survive.
while IFS= read -r -d '' datFile; do
    rel="${datFile#"$SRC"/}"          # e.g. GEPCellsTowersEtaSK/foo.dat
    subdir="${rel%%/*}"               # e.g. GEPCellsTowersEtaSK
    [ "$subdir" = "$rel" ] && subdir="."   # a .dat sitting at the top level

    settings_for "$subdir"

    outDir="$DST/$(dirname "$rel")"
    base="$(basename "$rel" .dat)"
    outTxt="$outDir/$base.txt"
    outDat="$outDir/$base.dat"

    maxArg=()
    [ "$maxWords" -gt 0 ] && maxArg=(--max-jets "$maxWords")

    echo "[$subdir] ${bits}-bit, maxWords=${maxWords:-all}: $rel"
    if [ "$DRY_RUN" -eq 1 ]; then
        echo "  would write $outTxt"
        continue
    fi

    mkdir -p "$outDir"
    if python3 "$CONVERTER" "$datFile" --bits "$bits" "${maxArg[@]}" -o "$outTxt"; then
        cp -f "$datFile" "$outDat"
        nConverted=$((nConverted + 1))
    else
        echo "  [error] conversion failed for $rel" >&2
        nFailed=$((nFailed + 1))
    fi
done < <(find "$SRC" -name '*.dat' -type f -print0 | sort -z)

echo
if [ "$DRY_RUN" -eq 1 ]; then
    echo "Dry run: nothing written. Destination would be $DST"
else
    echo "Converted $nConverted file(s) into $DST"
    [ "$nFailed" -gt 0 ] && echo "$nFailed file(s) FAILED -- see errors above" >&2
fi
exit $(( nFailed > 0 ? 1 : 0 ))
