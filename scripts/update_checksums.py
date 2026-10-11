#!/usr/bin/env python3
"""Refresh checksums after intentional source or data changes."""
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DIRECTORIES = ('src', 'drivers', 'examples', 'wrappers', 'recipes', 'patches',
               'third_party', 'provenance', 'data', 'scripts')


def main():
    files = [ROOT/'LICENSE', ROOT/'THIRD_PARTY_NOTICES.md']
    for directory in DIRECTORIES:
        files.extend(path for path in (ROOT/directory).rglob('*')
                     if path.is_file() and '__pycache__' not in path.parts
                     and path.suffix != '.pyc')
    lines = [hashlib.sha256(path.read_bytes()).hexdigest()+'  '+path.relative_to(ROOT).as_posix()
             for path in sorted(files)]
    (ROOT/'SHA256SUMS').write_text('\n'.join(lines)+'\n', encoding='utf-8', newline='\n')
    print('Updated SHA256SUMS for '+str(len(files))+' source, license, and data files.')


if __name__ == '__main__':
    main()
