"""First-mention monotonicity checker for NAR_MANUSCRIPT_REVISED.md.

Scans the body (before SUPPLEMENTARY MATERIAL) for Table N mentions (N=1..12),
expanding ranges like "Tables 1-11" and lists like "Tables 2 and 5".
Excludes "Supplementary Table" mentions. A table's own caption line
("**Table N.**" or "Table N." at start) counts as a normal occurrence.
Prints each table's first-mention line and whether the order is strictly 1..12.
"""
import re
import sys

PATH = r"D:/系统发育树项目/Fusang/Fusang-main/NAR_MANUSCRIPT_REVISED.md"


def main():
    lines = open(PATH, encoding="utf-8").read().split("\n")
    supp_idx = next(i for i, l in enumerate(lines)
                    if l.strip().upper().startswith("## SUPPLEMENTARY MATERIAL"))
    first = {}
    pat_single = re.compile(r"(?<!Supplementary )(?<!Supp )Table (\d+)")
    pat_range = re.compile(r"(?<!Supplementary )(?<!Supp )Tables (\d+)[\u2013-](\d+)")
    pat_list = re.compile(r"(?<!Supplementary )(?<!Supp )Tables (\d+(?: and \d+)+)")
    for ln, line in enumerate(lines[:supp_idx], start=1):
        for m in pat_range.finditer(line):
            a, b = int(m.group(1)), int(m.group(2))
            for n in (a, b):
                if 1 <= n <= 12 and n not in first:
                    first[n] = ln
        for m in pat_list.finditer(line):
            for tok in re.findall(r"\d+", m.group(1)):
                n = int(tok)
                if 1 <= n <= 12 and n not in first:
                    first[n] = ln
        line2 = pat_range.sub("", line)
        line2 = pat_list.sub("", line2)
        for m in pat_single.finditer(line2):
            n = int(m.group(1))
            if 1 <= n <= 12 and n not in first:
                first[n] = ln
    order = [n for n, _ in sorted(first.items(), key=lambda kv: kv[1])]
    for n in sorted(first):
        print(f"Table {n}: first mention at line {first[n]}")
    print("order:", order)
    print("monotonic 1..12:", order == list(range(1, 13)))
    return 0 if order == list(range(1, 13)) else 1


if __name__ == "__main__":
    sys.exit(main())
