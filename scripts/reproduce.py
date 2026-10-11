#!/usr/bin/env python3
"""Recompute the released feature-level evidence with Python's standard library."""
import argparse
from collections import defaultdict
import csv
import hashlib
import io
import json
import math
import os
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data'


def read(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))


def moments(values):
    n = len(values)
    if n < 2:
        return n, math.nan, math.nan
    mean = math.fsum(values) / n
    variance = math.fsum((x-mean)**2 for x in values) / (n-1)
    return n, mean, variance


def welch(fixed, random):
    nf, mf, vf = moments(fixed)
    nr, mr, vr = moments(random)
    if nf < 2 or nr < 2:
        return math.nan
    se2 = vf/nf + vr/nr
    return (mf-mr)/math.sqrt(se2) if se2 > 0 else math.nan


def cp_test(records, x='L0', y='L1'):
    mx = math.fsum(float(r[x]) for r in records)/len(records)
    my = math.fsum(float(r[y]) for r in records)/len(records)
    groups = defaultdict(list)
    for r in records:
        groups[r['group']].append((float(r[x])-mx)*(float(r[y])-my))
    return welch(groups['fixed'], groups['random'])


def grouped(records, key):
    result = defaultdict(list)
    for r in records:
        result[r[key]].append(r)
    return result


def feature_test(records, name):
    groups = defaultdict(list)
    for r in records:
        value = float(r[name])
        if math.isfinite(value):
            groups[r['group']].append(value)
    return welch(groups['fixed'], groups['random'])


def write(path, fields, rows):
    with path.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fields, lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def check_hashes():
    manifest = ROOT / 'SHA256SUMS'
    if not manifest.exists():
        raise RuntimeError('SHA256SUMS is missing')
    checked = 0
    for line in manifest.read_text(encoding='utf-8').splitlines():
        digest, name = line.split('  ', 1)
        actual = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        if actual != digest:
            raise RuntimeError('Checksum mismatch: '+name)
        checked += 1
    return checked


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-core', action='store_true', help='compile the extracted C reducers and compare against two-pass Python')
    parser.add_argument('--cc', default=os.environ.get('CC','gcc'), help='GCC/Clang executable; defaults to CC or gcc')
    args = parser.parse_args()
    files_checked = check_hashes()
    out = ROOT/'results'
    out.mkdir(exist_ok=True)
    protocol = json.loads((DATA/'protocol.json').read_text())
    metadata = {r['campaign']:r for r in read(DATA/'campaigns.csv')}
    outputs = grouped(read(DATA/'masked_outputs.csv'), 'campaign')
    calculated = []
    for campaign, records in outputs.items():
        assert len(records) == 5000 and len({int(r['trace_id']) for r in records}) == 5000
        assert sum(r['group']=='fixed' for r in records) == 2500
        for r in records:
            assert r['group'] == ('fixed' if int(r['trace_id'])%2==0 else 'random')
        calculated.append({'campaign':campaign,'n_fixed':2500,'n_random':2500,
                           't_hw0':feature_test(records,'L0'),'t_hw1':feature_test(records,'L1'),
                           't_cp':cp_test(records)})
    assert len(calculated) == 15
    write(out/'masked_statistics.csv', list(calculated[0]), calculated)

    profiles = grouped(read(DATA/'sec_add_entry_profiles.csv'), 'campaign')
    hd, figure = {}, []
    for campaign, records in profiles.items():
        assert len(records)==5000
        channel_tests = []
        for channel in ('g0','g1'):
            values = defaultdict(list)
            for r in records:
                values[r['group']].append(math.fsum(float(r[f'{channel}_e{i}']) for i in range(1,6))/5)
            channel_tests.append(welch(values['fixed'],values['random']))
        hd[campaign] = max(abs(t) for t in channel_tests)
        for entry in range(1,6):
            groups = defaultdict(list)
            for r in records:
                groups[r['group']].append(float(r[f'g0_e{entry}']))
            figure.append({'campaign':campaign,'entry':entry,
                           'delta_mean_hd_bits':moments(groups['fixed'])[1]-moments(groups['random'])[1]})
    write(out/'figure2_points.csv', list(figure[0]), figure)
    for campaign, records in grouped(read(DATA/'a2b_internal_hd.csv'),'campaign').items():
        assert len(records)==5000
        hd[campaign] = abs(feature_test(records,'raw0_hd'))
    for campaign, records in outputs.items():
        if metadata[campaign]['module']=='sec_and':
            assert all(float(r['internal_hd'])==0 for r in records)
            hd[campaign] = math.nan

    table2 = []
    for module in ('sec_and','sec_add','a2b_core'):
        selected = [r for r in calculated if metadata[r['campaign']]['module']==module and
                    metadata[r['campaign']]['condition'] in ('default','allones')]
        hw = [max(abs(r['t_hw0']),abs(r['t_hw1'])) for r in selected]
        cp = [abs(r['t_cp']) for r in selected]
        internal = [hd[r['campaign']] for r in selected]
        assert max(hw)<4.5
        assert max(cp)>=4.5 or any(t>=4.5 for t in internal)
        table2.append({'module':module,'HW_min':min(hw),'HW_max':max(hw),'CP_min':min(cp),'CP_max':max(cp),
                       'HD_min':min(internal),'HD_max':max(internal)})
    write(out/'table2.csv',list(table2[0]),table2)
    assert all(abs(r['t_cp'])>=4.5 for r in calculated if metadata[r['campaign']]['module']=='sec_add')

    table1 = []
    for kernel, records in grouped(read(DATA/'pqc_paired_comparisons.csv'),'kernel').items():
        coverage = protocol['pqc_coverage'][kernel]
        assert len(records) == coverage['positions']
        delta = max(float(r['reported_abs_delta']) for r in records)
        assert delta <= 5e-11
        table1.append({'kernel':kernel,'positions':len(records),'archived_max_abs_delta':delta,
                       'dedicated_traces':coverage['dedicated_traces']})
    write(out/'table1.csv',list(table1[0]),table1)
    assert all(float(r['paired_abs_delta_t'])<=2.37e-12 for r in read(DATA/'sec_add_paired_checks.csv'))

    intervention = defaultdict(list)
    for r in read(DATA/'intervention_traces.csv'):
        intervention[(r['variant'],r['seed'])].append(r)
    assert len(intervention)==6
    interface_tests = []
    archived = read(DATA/'intervention_statistics.csv')
    features = list(archived_row['feature'] for archived_row in archived[:9])
    for (variant,seed), records in intervention.items():
        assert len(records)==5000
        for r in records:
            assert int(r['cycles'])==26
            assert r['identity_ok']==r['functional_ok']==r['protocol_ok']=='1'
            x,y,s0,s1 = (int(r[k],16) for k in ('x','y','share0','share1'))
            assert (s0^s1)==(x&y)
            assert float(r['return_bus_valid_endpoint'])==(s0^s1).bit_count()
            assert float(r['share0_hw'])==s0.bit_count() and float(r['share1_hw'])==s1.bit_count()
            assert float(r['return_bus_adjacent_sum'])==s0.bit_count()+s1.bit_count()
            if variant=='baseline':
                assert float(r['result_reg_between_returns'])==(x&y).bit_count()
            else:
                assert float(r['result_reg_between_returns'])==0
                assert float(r['share_return0_write'])==s0.bit_count()
                assert float(r['share_return1_write'])==s1.bit_count()
        for feature in features+['output_share_cp']:
            t = cp_test(records,'share0_hw','share1_hw') if feature=='output_share_cp' else feature_test(records,feature)
            interface_tests.append({'variant':variant,'seed':seed,'feature':feature,'t_recomputed':t})
    # The archived CP feature label is retained instead of inferred from prose.
    cp_label = next(r['feature'] for r in archived if r['feature'] not in features)
    for r in interface_tests:
        label = cp_label if r['feature']=='output_share_cp' else r['feature']
        reference = next(a for a in archived if (a['variant'],a['seed'],a['feature'])==(r['variant'],r['seed'],label))
        expected = float(reference['t_offline'])
        actual = r['t_recomputed']
        assert math.isfinite(actual)==math.isfinite(expected)
        if math.isfinite(actual):
            assert abs(actual-expected)<1e-9, (r,expected)
    for seed in ('0','1','2'):
        a,b = intervention[('baseline',seed)],intervention[('split_storage',seed)]
        for x,y in zip(a,b):
            for key in ('trace_id','group','input_hash','output_hash','protocol_hash','cycles','x','y','share0','share1'):
                assert x[key]==y[key], (seed,key)
    write(out/'intervention_recomputed.csv',list(interface_tests[0]),interface_tests)

    core_delta = None
    if args.check_core:
        cc = shutil.which(args.cc)
        if cc is None:
            raise RuntimeError('Compiler not found: '+args.cc)
        compiler_env = os.environ.copy()
        compiler_env['PATH'] = str(Path(cc).resolve().parent) + os.pathsep + compiler_env.get('PATH','')
        build = ROOT/'build'
        build.mkdir(exist_ok=True)
        exe = build/('replay_features.exe' if os.name=='nt' else 'replay_features')
        subprocess.run([cc,'-std=c11','-O2',str(ROOT/'examples/replay_features.c'),'-o',str(exe),'-lm'],check=True,env=compiler_env)
        text = subprocess.check_output([str(exe),str(DATA/'masked_outputs.csv')],text=True,env=compiler_env)
        actual = {r['campaign']:r for r in csv.DictReader(io.StringIO(text))}
        assert set(actual)==set(outputs)
        core_delta = 0.0
        for r in calculated:
            for key in ('t_hw0','t_hw1','t_cp'):
                error = abs(float(actual[r['campaign']][key])-r[key])
                assert error<1e-9, (r['campaign'],key,error)
                core_delta = max(core_delta,error)
        (out/'core_replay.csv').write_text(text,encoding='utf-8',newline='\n')

    report = {'status':'PASS','checksummed_files':files_checked,'output_campaigns':15,
              'masked_output_rows':sum(map(len,outputs.values())),'sec_add_profile_rows':sum(map(len,profiles.values())),
              'interface_campaigns':6,'interface_trace_rows':30000,'modules_with_indications_beyond_output_HW':3,
              'core_replay_max_abs_delta':core_delta,
              'scope':'feature-level recomputation and archived PQC residual checks; RTL execution is a separate entry point'}
    (out/'verification.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))
    print('Generated table1.csv, table2.csv, figure2_points.csv, and intervention_recomputed.csv in results/.')


if __name__=='__main__':
    main()
