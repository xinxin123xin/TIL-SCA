#!/usr/bin/env python3
"""Build one actual RTL experiment from a separately obtained WPI source tree."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def apply_patch(target, patch):
    """Apply our single-file unified patches with exact context, without Git."""
    original = target.read_text(encoding='utf-8').splitlines(keepends=True)
    lines = patch.read_text(encoding='utf-8').splitlines(keepends=True)
    result, cursor, i = [], 0, 2
    while i < len(lines):
        match = re.match(r'@@ -(\d+)(?:,\d+)? \+\d+(?:,\d+)? @@',lines[i])
        if not match:
            raise RuntimeError('Invalid patch hunk: '+str(patch))
        start = int(match.group(1))-1
        result.extend(original[cursor:start])
        cursor = start
        i += 1
        while i < len(lines) and not lines[i].startswith('@@ '):
            marker, text = lines[i][0], lines[i][1:]
            if marker in (' ','-'):
                if cursor >= len(original) or original[cursor] != text:
                    raise RuntimeError('Patch context mismatch: '+str(target))
                cursor += 1
            if marker in (' ','+'):
                result.append(text)
            elif marker not in ('-','\\'):
                raise RuntimeError('Unsupported patch line: '+lines[i])
            i += 1
    result.extend(original[cursor:])
    target.write_text(''.join(result),encoding='utf-8',newline='\n')


def run(command, work, log, env=None):
    with (work/log).open('w',encoding='utf-8') as stream:
        completed = subprocess.run(command,cwd=work,stdout=stream,stderr=subprocess.STDOUT,env=env)
    if completed.returncode:
        print((work/log).read_text(encoding='utf-8',errors='replace')[-4000:])
        raise RuntimeError('Build failed; see '+str(work/log))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',required=True,choices=['sec_and','sec_add','sec_add_entries','a2b','a2b_entries','interface_baseline','interface_split'])
    parser.add_argument('--upstream',required=True,type=Path,help='WPI repository/archive root containing rtl/obi and rtl/ise')
    parser.add_argument('--yosys',default='yosys')
    parser.add_argument('--cc',default=os.environ.get('CC','gcc'),help='GCC/Clang executable; defaults to CC or gcc')
    parser.add_argument('--condition',choices=['allones','zero','mixed'],default='allones',help='sec_add_entries fixed-operand condition')
    args = parser.parse_args()
    tools = {name:shutil.which(value) for name,value in [('yosys',args.yosys),('cc',args.cc)]}
    if not all(tools.values()):
        raise RuntimeError('Yosys and a GCC/Clang C compiler must be on PATH, or specified explicitly')
    compiler_env = os.environ.copy()
    compiler_bin = str(Path(tools['cc']).resolve().parent)
    compiler_env['PATH'] = compiler_bin + os.pathsep + compiler_env.get('PATH','')
    interface = args.case.startswith('interface_')
    case = 'interface' if interface else args.case.replace('_entries','')
    work = ROOT/'build'/(args.case+'_'+args.condition if args.case=='sec_add_entries' else args.case)
    work.mkdir(parents=True,exist_ok=True)
    target_dir = work/'rtl' if interface else work
    target_dir.mkdir(exist_ok=True)
    specs = json.loads((ROOT/'provenance/rtl_sources.json').read_text(encoding='utf-8'))[case]
    for spec in specs:
        upstream = args.upstream/spec['upstream_path']
        text = upstream.read_text(encoding='utf-8')
        if hashlib.sha256(text.encode()).hexdigest()!=spec['base_text_sha256']:
            raise RuntimeError('Upstream file differs from the evaluated snapshot: '+spec['upstream_path'])
        target = target_dir/spec['filename']
        target.write_text(text,encoding='utf-8',newline='\n')
        if spec['patch']:
            apply_patch(target,ROOT/spec['patch'])
        if hashlib.sha256(target.read_text(encoding='utf-8').encode()).hexdigest()!=spec['evaluated_text_sha256']:
            raise RuntimeError('Adapted RTL checksum mismatch: '+spec['filename'])
    if args.case=='interface_split':
        apply_patch(target_dir/'sca_pqc_coprocessor_top.sv',ROOT/'patches/split_return_storage.patch')
    if not interface:
        for wrapper in (ROOT/'wrappers'/case).glob('*.sv'):
            shutil.copyfile(wrapper,work/wrapper.name)
    shutil.copyfile(ROOT/'third_party/wpi/tsim.h',work/'tsim.h')
    shutil.copyfile(ROOT/'recipes'/f'{case}.ys',work/'generate.ys')
    run([tools['yosys'],'-s','generate.ys'],work,'yosys.log')
    flat = {'sec_and':'dom_and_flat.json','sec_add':'sec_add_flat.json','a2b':'a2b_flat.json','interface':'coproc_flat.json'}[case]
    run([sys.executable,str(ROOT/'third_party/wpi/gen_oim_fullchip.py'),'--json',flat,
         '--out-h','oim_bitlvl.h','--out-json','oim_bitlvl.json'],work,'oim.log')
    driver = {'sec_and':'sec_and_driver.c','sec_add':'sec_add_output_driver.c','sec_add_entries':'sec_add_entry_driver.c',
              'a2b':'a2b_output_driver.c','a2b_entries':'a2b_entry_driver.c',
              'interface_baseline':'intervention_driver.c','interface_split':'intervention_driver.c'}[args.case]
    defines = []
    if args.case=='sec_add_entries':
        a,b = {'allones':('0x1fff','0x1fff'),'zero':('0','0'),'mixed':('0x0a53','0x136c')}[args.condition]
        defines = ['-DFIXED_A='+a,'-DFIXED_B='+b]
    exe = 'driver.exe' if sys.platform=='win32' else 'driver'
    run([tools['cc'],'-std=c11','-O2','-I.',*defines,str(ROOT/'drivers'/driver),'-o',exe,'-lm'],work,'compile.log',env=compiler_env)
    print('Compiler: '+str(tools['cc'])+' (bin directory prepended to build PATH)')
    print('Built '+str(work/exe))
    print('Run from this directory so generated CSV files stay in build/:')
    print('  cd "'+str(work)+'"')
    print('  '+('.\\' if sys.platform=='win32' else './')+exe+' 10'+(' --campaign-seed 0' if case=='sec_and' else ' 0'))
    print('10 samples per class is a functional/numerical smoke test; paper campaigns use 2500 per class.')


if __name__=='__main__':
    main()
