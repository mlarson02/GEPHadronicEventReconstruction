#!/usr/bin/env python3
"""Convert jet .dat files to space-separated hex word .txt format.

Each output line corresponds to one event. Hex values are extracted from the
trailing 0xNNN field of each jet entry and zero-padded to the specified bit
width (64 or 128 bits).

Usage:
    python3 convert_dat.py input_data_seed.dat --bits 128
    python3 convert_dat.py input_data.dat --bits 64 -o input01_advanced.txt
"""

import argparse
import os
import re
import sys


def convert_dat(input_file: str, bits: int = 128, max_jets: int | None = None,
                output_file: str | None = None) -> None:
    hex_chars = bits // 4  # 32 for 128-bit, 16 for 64-bit
    trailing_hex = re.compile(r'0x([0-9a-fA-F]+)\s*$')

    events: list[list[str]] = []
    current: list[str] | None = None

    with open(input_file, 'r') as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('Event :'):
                if current is not None:
                    events.append(current)
                current = []
            elif current is not None:
                m = trailing_hex.search(line)
                if m:
                    current.append(m.group(1).zfill(hex_chars))

    if current is not None:
        events.append(current)

    if not events:
        print(f"No events found in {input_file}", file=sys.stderr)
        sys.exit(1)

    if output_file is None:
        base, _ = os.path.splitext(input_file)
        output_file = base + '.txt'

    with open(output_file, 'w') as f:
        for event in events:
            words = event[:max_jets] if max_jets is not None else event
            f.write(' '.join(words) + '\n')

    max_seen = max(len(e) for e in events)
    print(f"Wrote {len(events)} events (max {max_seen} words/line, {bits}-bit) -> {output_file}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Convert jet .dat file to space-separated hex .txt format')
    parser.add_argument('input', help='Input .dat file path')
    parser.add_argument('--bits', type=int, choices=[64, 128], default=128,
                        help='Bit width for zero-padding each word (default: 128)')
    parser.add_argument('--max-jets', type=int, default=None,
                        help='Truncate each event to at most N words (no zero-padding)')
    parser.add_argument('-o', '--output', default=None,
                        help='Output .txt file path (default: replaces .dat extension)')
    args = parser.parse_args()
    convert_dat(args.input, args.bits, args.max_jets, args.output)


if __name__ == '__main__':
    main()
