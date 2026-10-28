"""
regions.py — reading BED-like region files (motif/ChIP intersections, BED6, refFlat TSS).

Column layout is sniffed from the column count of the first data line:
  - 22 cols: peak+motif intersection -> chr,start,end,strand at 13,14,15,16
    (the 11-col motif layout, offset by an 11-col peak summary in front of it)
  - 11 cols: motif/ChIP intersection  -> chr,start,end,strand at 2,3,4,5
  -  6 cols: standard BED6            -> chr,start,end,strand at 0,1,2,5
  -  7 cols: refFlat TSS bed          -> chr,left,right,strand at 0,1,2,3
  -  4 cols: chrom,start,end,strand   -> 0,1,2,3
Any other column count raises instead of being guessed at.

Region centre is (start + end) // 2 everywhere (metaplots and logos), which
assumes the BED interval's midpoint lands on the motif centre. Fix the BED, not
the code, if it doesn't.
"""

import gzip
import os


def open_regions(path):
    """Open a region file as text, detecting gzip by magic bytes rather than the .gz suffix."""
    with open(path, 'rb') as fh:
        is_gzipped = fh.read(2) == b'\x1f\x8b'
    return gzip.open(path, 'rt') if is_gzipped else open(path)


def first_data_line(path):
    with open_regions(path) as fh:
        for line in fh:
            if line.strip() and not line.startswith('#'):
                return line
    return ''


def label_from_basename(regions_path):
    """Region file's basename with any .gz and then any single extension
    stripped, sanitised to filename-safe characters.
    'TSS_K562_RPKM.Q1.bed' -> 'TSS_K562_RPKM.Q1'."""
    base = os.path.basename(regions_path)
    if base.endswith('.gz'):
        base = base[:-3]
    stem = os.path.splitext(base)[0]
    safe = ''.join(c if (c.isalnum() or c in '._-') else '_' for c in stem)
    return safe or 'SIGNAL'


def region_layout(regions_path):
    """Returns (field_indices, analysis_label), field_indices being
    (chr_idx, start_idx, end_idx, strand_idx).

    11/22-col files are labelled by their motif-name column; everything else by
    the file's own basename stem, which is unique per file so outputs from
    different region files cannot overwrite each other."""
    fields = first_data_line(regions_path).rstrip('\r\n').split('\t')
    n = len(fields)
    if n in (11, 22):
        offset = n - 11
        label = fields[1 + offset] or 'MOTIF'
        return (2 + offset, 3 + offset, 4 + offset, 5 + offset), label
    label = label_from_basename(regions_path)
    if n == 6:
        # BED6: strand is column 5; column 3 is the name
        return (0, 1, 2, 5), label
    if n in (4, 7):
        # 7: refFlat TSS bed, '#chr left right strand geneName(s) geneID(s) transcript(s)'
        return (0, 1, 2, 3), label
    raise ValueError(
        f'{regions_path}: unrecognised column count {n}. Supported layouts are '
        '4 (chrom,start,end,strand), 6 (BED6), 7 (refFlat TSS: '
        'chr,left,right,strand), 11 (motif/ChIP intersection) and 22 '
        '(peak+motif intersection).')


def parse_region_line(line, field_indices):
    """(chrom, start, end, strand) strings from one data line. The line ending is
    stripped first: in the 4/6-col layouts strand is the last field and would
    otherwise come back as '+\\n'."""
    fields = line.rstrip('\r\n').split('\t')
    ci, si, ei, sti = field_indices
    return fields[ci], fields[si], fields[ei], fields[sti]


def iter_region_tuples(regions_path):
    """(chrom, start, end, strand) generator over every data line."""
    field_indices, _ = region_layout(regions_path)
    with open_regions(regions_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            yield parse_region_line(line, field_indices)
