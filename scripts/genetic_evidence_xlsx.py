"""Read exact XLSX cells with stdlib, keeping worksheet/cell coordinates."""
import posixpath
import xml.etree.ElementTree as ET
import zipfile
NS = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
class Workbook:
    def __init__(self, path):
        self.z = zipfile.ZipFile(path)
        self.strings = []
        if 'xl/sharedStrings.xml' in self.z.namelist():
            with self.z.open('xl/sharedStrings.xml') as f:
                for _, el in ET.iterparse(f, events=('end',)):
                    if el.tag == '{'+NS['m']+'}si':
                        self.strings.append(''.join(x.text or '' for x in el.iter('{'+NS['m']+'}t')))
                        el.clear()
        rels={x.attrib['Id']:x.attrib['Target'] for x in ET.fromstring(self.z.read('xl/_rels/workbook.xml.rels'))}
        self.sheets={x.attrib['name']:posixpath.normpath('xl/'+rels[x.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']]).lstrip('/') for x in ET.fromstring(self.z.read('xl/workbook.xml')).findall('m:sheets/m:sheet',NS)}
    def rows(self, sheet):
        with self.z.open(self.sheets[sheet]) as f:
            for _, el in ET.iterparse(f,events=('end',)):
                if el.tag != '{'+NS['m']+'}row': continue
                cells={}
                for c in el.findall('m:c',NS):
                    v=c.find('m:v',NS)
                    value=v.text if v is not None else ''
                    if c.attrib.get('t')=='s': value=self.strings[int(value)]
                    elif c.attrib.get('t')=='inlineStr': value=''.join(x.text or '' for x in c.findall('.//m:t',NS))
                    if value: cells[c.attrib['r']]=value
                yield int(el.attrib['r']),cells
                el.clear()
