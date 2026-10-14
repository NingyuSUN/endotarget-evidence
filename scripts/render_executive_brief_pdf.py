#!/usr/bin/env python3
"""Render the English brief as exactly two A4 pages without dropping source content."""
import argparse,hashlib,html,json,pathlib,re,datetime,importlib.metadata
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,PageBreak,Table,TableStyle
from pypdf import PdfReader
root=pathlib.Path(__file__).resolve().parent.parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--receipt', default='receipts/executive_brief_pdf.json')
parser.add_argument('--preview-dir', default='work/local-runtime/pdf-preview')
args=parser.parse_args()
source=root/"reports/executive_brief_en.md"
output=root/"reports/executive_brief_en.pdf"
font_dir=pathlib.Path("C:/Windows/Fonts")
if font_dir.exists():
    files=("arial.ttf","arialbd.ttf","ariali.ttf")
else:
    font_dir=pathlib.Path("/usr/share/fonts/truetype/dejavu")
    files=("DejaVuSans.ttf","DejaVuSans-Bold.ttf","DejaVuSans-Oblique.ttf")
for name,file in zip(("Brief","Brief-Bold","Brief-Italic"),files):
    pdfmetrics.registerFont(TTFont(name,str(font_dir/file)))
pdfmetrics.registerFontFamily("Brief",normal="Brief",bold="Brief-Bold",italic="Brief-Italic",boldItalic="Brief-Bold")
raw=source.read_text(encoding="utf-8")
blue=colors.HexColor("#17334A"); muted=colors.HexColor("#496372")
styles={
 "body":ParagraphStyle("body",fontName="Brief",fontSize=9.4,leading=12.6,textColor=blue,spaceAfter=5),
 "title":ParagraphStyle("title",fontName="Brief-Bold",fontSize=18,leading=21,textColor=blue,spaceAfter=11),
 "heading":ParagraphStyle("heading",fontName="Brief-Bold",fontSize=11,leading=14,textColor=blue,spaceBefore=7,spaceAfter=5),
 "table":ParagraphStyle("table",fontName="Brief",fontSize=8.4,leading=11,textColor=blue),
 "tablehead":ParagraphStyle("tablehead",fontName="Brief-Bold",fontSize=8.2,leading=10.7,textColor=blue),
}
def markup(text):
    result=html.escape(text)
    result=re.sub(r"\[([^\]]+)\]\(([^)]+)\)",lambda m:'<a href="'+m.group(2)+'" color="#176989">'+m.group(1)+'</a>',result)
    result=re.sub(r"\*\*(.+?)\*\*",r"<b>\1</b>",result)
    result=re.sub(r"\x60([^\x60]+)\x60",r"\1",result)
    return result
def paragraph(s,style="body"): return Paragraph(markup(s),styles[style])
lines=raw.splitlines(); story=[]; index=0
while index<len(lines):
    line=lines[index].strip()
    if not line:
        index+=1; continue
    if line.startswith("## Evidence integration"):
        story.append(PageBreak())
    if line.startswith("# "):
        story.append(paragraph(line[2:],"title")); index+=1; continue
    if line.startswith("## "):
        story.append(paragraph(line[3:],"heading")); index+=1; continue
    if line.startswith("|"):
        rows=[]
        while index<len(lines) and lines[index].strip().startswith("|"):
            cells=[x.strip() for x in lines[index].strip().strip("|").split("|")]
            if not all(re.fullmatch(r":?-+:?",x.replace(" ","")) for x in cells):
                rows.append(cells)
            index+=1
        width=A4[0]-84
        ratios=[.28,.25,.19,.14,.14] if len(rows[0])==5 else [.28,.25,.24,.23]
        values=[[paragraph(c,"tablehead" if n==0 else "table") for c in row] for n,row in enumerate(rows)]
        table=Table(values,colWidths=[width*r for r in ratios],hAlign="LEFT",repeatRows=1)
        table.setStyle(TableStyle([("BACKGROUND",(0,0),(-1,0),colors.HexColor("#EAF1F5")),
           ("VALIGN",(0,0),(-1,-1),"TOP"),("LEFTPADDING",(0,0),(-1,-1),5),
           ("RIGHTPADDING",(0,0),(-1,-1),5),("TOPPADDING",(0,0),(-1,-1),5),
           ("BOTTOMPADDING",(0,0),(-1,-1),5),
           ("LINEBELOW",(0,-1),(-1,-1),.5,colors.HexColor("#CBD8E1"))]))
        story.extend([table,Spacer(1,7)]); continue
    parts=[line]; index+=1
    while index<len(lines) and lines[index].strip() and not lines[index].lstrip().startswith(("#","|")):
        parts.append(lines[index].strip());index+=1
    story.append(paragraph(" ".join(parts)))
def footer(canvas,doc):
    canvas.setFont("Brief",8);canvas.setFillColor(muted)
    canvas.drawString(42,25,"EndoTarget Evidence  |  Research decision brief")
    canvas.drawRightString(A4[0]-42,25,str(doc.page)+" / 2")
doc=SimpleDocTemplate(str(output),pagesize=A4,rightMargin=42,leftMargin=42,
  topMargin=35,bottomMargin=42,title="EndoTarget Evidence — research decision brief",
  author="EndoTarget Evidence project",pageCompression=1)
doc.build(story,onFirstPage=footer,onLaterPages=footer)
reader=PdfReader(str(output)); pages=[page.extract_text() for page in reader.pages]
if len(pages)!=2: raise SystemExit(f"Expected exactly 2 pages, got {len(pages)}; adjust layout, never omit text")
for value in ["150,788","52,744","487","13","SKIDA1","PARP9","unavailable","500"]:
    if value not in "\n".join(pages): raise SystemExit("Missing PDF evidence: "+value)
# Verify all plain paragraph/table text tokens survived extraction; allow PDF whitespace.
plain=re.sub(r"\[([^\]]+)\]\([^)]+\)",r"\1",raw)
plain=plain.replace("**","").replace("#","").replace("|"," ")
tokens=[x for x in re.findall(r"[A-Za-z0-9][A-Za-z0-9,_]+",plain) if not re.fullmatch(r"[-:]+",x)]
rendered=" ".join(pages)
missing=sorted({x for x in tokens if x not in rendered})
if missing: raise SystemExit("PDF text omissions: "+repr(missing))
receipt={"created_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
 "source":source.relative_to(root).as_posix(),"source_sha256":hashlib.sha256(source.read_bytes()).hexdigest(),
 "output":output.relative_to(root).as_posix(),"output_sha256":hashlib.sha256(output.read_bytes()).hexdigest(),
 "script_sha256":hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),
 "pages":len(pages),"text_coverage":"All source alphanumeric tokens and key scientific values checked",
 "new_scientific_analysis":False,"visual_review":"PENDING",
 "renderer_versions":{name:importlib.metadata.version(name) for name in ("reportlab","pypdf")},
 "font_sha256":{file:hashlib.sha256((font_dir/file).read_bytes()).hexdigest() for file in files}}
preview=root/args.preview_dir;preview.mkdir(parents=True,exist_ok=True)
try:
    import pypdfium2
    pdf=pypdfium2.PdfDocument(str(output))
    for i,page in enumerate(pdf):
        page.render(scale=1.5).to_pil().save(preview/f"brief-{i+1}.png")
    receipt["preview_rendered"]=True
except ImportError:
    receipt["preview_rendered"]=False
(root/args.receipt).write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
print(json.dumps(receipt))
