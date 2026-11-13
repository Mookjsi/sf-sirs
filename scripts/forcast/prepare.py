"""Extract RSV weekly counts from local annual and IDWR reports."""
import argparse, csv, hashlib, re, zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

NS = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
JAPANESE_PREFECTURES = '''北海道 青森県 岩手県 宮城県 秋田県 山形県 福島県 茨城県 栃木県 群馬県 埼玉県 千葉県 東京都 神奈川県 新潟県 富山県 石川県 福井県 山梨県 長野県 岐阜県 静岡県 愛知県 三重県 滋賀県 京都府 大阪府 兵庫県 奈良県 和歌山県 鳥取県 島根県 岡山県 広島県 山口県 徳島県 香川県 愛媛県 高知県 福岡県 佐賀県 長崎県 熊本県 大分県 宮崎県 鹿児島県 沖縄県'''.split()
ENGLISH_PREFECTURES = '''Hokkaido Aomori Iwate Miyagi Akita Yamagata Fukushima Ibaraki Tochigi Gunma Saitama Chiba Tokyo Kanagawa Niigata Toyama Ishikawa Fukui Yamanashi Nagano Gifu Shizuoka Aichi Mie Shiga Kyoto Osaka Hyogo Nara Wakayama Tottori Shimane Okayama Hiroshima Yamaguchi Tokushima Kagawa Ehime Kochi Fukuoka Saga Nagasaki Kumamoto Oita Miyazaki Kagoshima Okinawa'''.split()
PREFECTURES = dict(zip(JAPANESE_PREFECTURES, ENGLISH_PREFECTURES))

def cells(data, strings):
    root = ET.fromstring(data)
    result = {}
    for cell in root.findall('.//m:sheetData/m:row/m:c', NS):
        value = cell.find('m:v', NS)
        inline = cell.find('m:is', NS)
        text = value.text if value is not None else ''.join(inline.itertext()) if inline is not None else ''
        if cell.get('t') == 's': text = strings[int(text)]
        result[cell.get('r')] = text
    return result

def write(path, columns, rows):
    with open(path, 'w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader(); writer.writerows(rows)

def weekly_idwr(path, population):
    with open(path, encoding='cp932', newline='') as stream:
        table = list(csv.reader(stream))
    header = table[1][0]
    match = re.match(r'(20\d\d)年(\d\d)週', header)
    if not match:
        raise ValueError(f'{path}: expected year/week in second row, got {header!r}')
    year, last_week = map(int, match.groups())
    starts = [i for i, row in enumerate(table) if row and row[0] == 'ＲＳウイルス感染症']
    if len(starts) != 1:
        raise ValueError(f'{path}: expected one RSV disease block, found {len(starts)}')
    start = starts[0]
    weeks = [int(table[start+1][col].replace('週', '')) for col in range(3, len(table[start+1]), 2)
             if table[start+1][col]]
    if weeks != list(range(1, last_week+1)):
        missing = sorted(set(range(1, last_week+1))-set(weeks))
        raise ValueError(f'{path}: missing weekly columns {missing}; supply files ' +
                         ', '.join(f'{year}-{w:02d}-teiten-tougai.csv' for w in missing))
    if table[start+2][1:5] != ['報告', '定当', '報告', '定当']:
        raise ValueError(f'{path}: expected cumulative and weekly report/rate column pairs')
    names = [table[start+4+i][0] for i in range(47)]
    if names != JAPANESE_PREFECTURES:
        raise ValueError(f'{path}: expected 47 prefectures in official order; got {names!r}')
    if set(ENGLISH_PREFECTURES) != set(population):
        raise ValueError('Population mapping does not match the 47 official prefectures')
    rows = []
    for i, name in enumerate(names):
        raw = table[start+4+i]
        values = []
        for week in weeks:
            col = 3+2*(week-1)
            if table[start+1][col] != f'{week:02d}週' or table[start+2][col] != '報告':
                raise ValueError(f'{path}: weekly count header mismatch at W{week:02d}')
            value = raw[col].strip() if col < len(raw) else ''
            count = int(value) if re.fullmatch(r'\d+', value) else None
            values.append(count)
            rows.append(dict(year=year, week=week, prefecture=PREFECTURES[name],
                             cases='' if count is None else count, island=population[PREFECTURES[name]],
                             source_version='JIHS IDWR weekly CSV', source_path=str(path.resolve()),
                             missing_reason='' if count is not None else f'blank/non-numeric report: {value!r}'))
        cumulative = raw[1].strip()
        if all(v is not None for v in values) and (not cumulative.isdigit() or sum(values) != int(cumulative)):
            raise ValueError(f'{path}: RSV cumulative count disagrees with weekly counts for {name}')
    for week in weeks:
        reports = [int(r['cases']) for r in rows if r['week'] == week and r['cases'] != '']
        national = table[start+3][3+2*(week-1)].strip()
        if len(reports) == 47 and (not national.isdigit() or sum(reports) != int(national)):
            raise ValueError(f'{path}: RSV national total disagrees at W{week:02d}')
    return year, last_week, rows

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--annual', default='data/forcast/Syu_08_1.xlsx')
    p.add_argument('--population', default='data/processed/original-island-v1-traceable/japan_pop_prefecture.csv')
    p.add_argument('--cases', default='data/processed/original-island-v1-traceable/japan_rsv_prefecture.csv')
    p.add_argument('--out', required=True)
    p.add_argument('--weekly', action='append', type=Path,
                   help='Local year-to-date IDWR weekly CSV; repeat for each year')
    a = p.parse_args(); out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    for path in (a.annual, a.population, a.cases):
        if not Path(path).is_file(): raise SystemExit(f'Missing input: {path}')
    with open(a.population, newline='') as stream: population = list(csv.DictReader(stream))
    honshu = {r['prefecture'] for r in population if r['island'] == 'Honshu'}
    assert len(population) == 47 and len(honshu) == 34
    population_map = {r['prefecture']: r['island'] for r in population}
    with zipfile.ZipFile(a.annual) as book:
        strings = []
        if 'xl/sharedStrings.xml' in book.namelist():
            strings = [''.join(x.itertext()) for x in ET.fromstring(book.read('xl/sharedStrings.xml'))]
        sheets = sorted((n for n in book.namelist() if re.fullmatch(r'xl/worksheets/sheet\d+\.xml', n)), key=lambda n:int(re.search(r'sheet(\d+)',n).group(1)))
        assert len(sheets) >= 53, f'Expected annual total plus 52 weekly sheets, found {len(sheets)}'
        annual = cells(book.read(sheets[0]), strings)
        assert '2024' in annual.get('B2','') and 'Respiratory syncytial virus' in annual.get('I5','') and 'total No.' in annual.get('I6','')
        rows = []; national = {}
        for week in range(1,53):
            sheet = cells(book.read(sheets[week]), strings)
            assert '2024' in sheet.get('B2','') and re.search(r'week\s+'+str(week)+r'\)',sheet.get('H3',''), re.I), (sheets[week],sheet.get('H3'))
            assert 'Respiratory syncytial virus' in sheet.get('I5','') and 'total No.' in sheet.get('I6','')
            seen = set(); total = 0
            for row in range(8,55):
                label = sheet.get(f'B{row}','')
                match = re.search(r'\(([^()]*)\)',label)
                assert match, (week,row,label)
                prefecture = match.group(1).strip()
                assert prefecture in {x['prefecture'] for x in population} and prefecture not in seen, (week,prefecture)
                seen.add(prefecture)
                raw = sheet.get(f'I{row}','')
                assert re.fullmatch(r'\d+',raw), f'Non-numeric or missing report count: week {week}, {prefecture}, {raw!r}'
                count = int(raw); total += count
                rows.append(dict(year=2024,week=week,prefecture=prefecture,cases=count,island='Honshu' if prefecture in honshu else 'Other',source_version='JIHS finalized NESID annual table 8-1'))
            assert len(seen)==47 and total==int(sheet['I7']), (week,total,sheet.get('I7'))
            national[week]=total
        assert sum(national.values())==int(annual['I7'])
    for row in rows:
        row['source_path'] = str(Path(a.annual).resolve())
        row['missing_reason'] = ''
    weekly_paths = a.weekly if a.weekly is not None else sorted(Path('data/forcast').glob('*-teiten-tougai.csv'))
    by_year = {}
    for path in weekly_paths:
        if not path.is_file(): raise SystemExit(f'Missing weekly IDWR CSV: {path}')
        year, last_week, new_rows = weekly_idwr(path, population_map)
        if year in by_year:
            raise SystemExit(f'Duplicate year {year} in weekly IDWR inputs; supply one latest file per year')
        by_year[year] = (last_week, new_rows)
        rows.extend(new_rows)
    keys = [(r['year'], r['week'], r['prefecture']) for r in rows]
    if len(keys) != len(set(keys)):
        raise SystemExit('Duplicate year/week/prefecture observation key')
    write(out/'honshu_external_cases.csv',
          ['year','week','prefecture','cases','island','source_version','source_path','missing_reason'], rows)
    missing_weeks = []
    for year, (last_week, year_rows) in by_year.items():
        for week in range(1, last_week+1):
            present = {r['prefecture'] for r in year_rows if r['week'] == week and r['cases'] != ''}
            if len(present & honshu) != 34:
                missing_weeks.append(dict(year=year, week=week, reason='missing Honshu prefecture report',
                    missing_prefectures=';'.join(sorted(honshu-present))))
    write(out/'observation_exclusions.csv',['year','week','reason','missing_prefectures'],missing_weeks)
    with open(a.cases,newline='') as stream: history=list(csv.DictReader(stream))
    prior=[r for r in history if int(r['year'])==2023 and r['prefecture'] in honshu]
    assert len(prior)==34*52 and len({(r['week'],r['prefecture']) for r in prior})==len(prior)
    weekly={w:sum(int(r['cases']) for r in prior if int(r['week'])==w) for w in range(1,53)}
    assert weekly[52]==120 and sum(weekly.values())==102497 and max(weekly,key=weekly.get)==27 and weekly[27]==7785
    national_prior=sum(int(r['cases']) for r in history if int(r['year'])==2023 and int(r['week'])==52)
    assert national_prior==189
    write(out/'naive_2023.csv',['week','cases'],[dict(week=w,cases=weekly[w]) for w in weekly])
    write(out/'overlap_audit.csv',['status','validated_weeks','note'],[dict(status='not_performed',validated_weeks=0,note='No optional 2023 IDWR weekly CSV supplied; competition-file checks passed')])
    manifest=[]
    for path,kind in ((a.annual,'finalized_annual_2024'),(a.population,'population_mapping'),(a.cases,'competition_cases_2013_2023')) + tuple((str(path),'IDWR_weekly_CSV') for path in weekly_paths):
        manifest.append(dict(path=str(Path(path).resolve()),kind=kind,sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest()))
    write(out/'source_manifest.csv',['path','kind','sha256'],manifest)
    print(f'Validated 2024 annual plus {[(year, item[0]) for year,item in by_year.items()]}; incomplete weeks {len(missing_weeks)}; 2024 national sum {sum(national.values())}; W52 {national[52]}')
if __name__=='__main__': main()
