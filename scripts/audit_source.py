"""Inventory hand-written production MoonBit and basic documentation integrity.
Physical lines include normal in-file documentation; code-only is also reported.
No tests/examples/generated/vendor paths are counted towards the production gate.
"""
import json
import re
import subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def audit():
    files=[]
    for path in sorted(ROOT.glob('*.mbt')):
        if path.name.endswith(('_test.mbt','_wbtest.mbt')):
            continue
        text=path.read_text(encoding='utf-8-sig')
        lines=text.splitlines()
        assert not re.search(r'\b(?:TODO|FIXME|abort|panic)\s*\(', text), path.name
        files.append({'file':path.name,'physical_lines':len(lines),'code_lines':sum(bool(x.strip()) and not x.lstrip().startswith('//') for x in lines)})
    total=sum(x['physical_lines'] for x in files)
    assert total>=1500, f'production physical lines {total}<1500'
    tests=sum(len(re.findall(r'(?m)^test\s+"',p.read_text(encoding='utf-8-sig'))) for p in ROOT.glob('*test.mbt'))
    for path in [ROOT/'README.md',ROOT/'docs'/'architecture.md',ROOT/'docs'/'validation'/'independent-oracle.md']:
        for target in re.findall(r'\]\(([^)]+)\)',path.read_text(encoding='utf-8-sig')):
            if not target.startswith(('https://','http://','#')):
                assert (path.parent/target.split('#')[0]).is_file(), f'broken link {path}: {target}'
    commits=subprocess.check_output(['git','rev-list','--count','HEAD'],cwd=ROOT,text=True).strip()
    return {'production_files':files,'production_physical_lines':total,'production_code_lines':sum(x['code_lines'] for x in files),'test_cases':tests,'raw_commits':int(commits),'commit_quality':'requires human semantic audit; raw count alone is not acceptance'}
if __name__=='__main__':
    print(json.dumps(audit(),ensure_ascii=False,indent=2))
