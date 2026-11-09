from pathlib import Path
import zipfile, math, copy
import xml.etree.ElementTree as E
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / '수정본_20261103'
SRC = BASE / 'UAM_Ring_Candidate_Selection_Process.pptx'
OUT = BASE / 'UAM_Landing_Area_Vertical_Penalty.pptx'
TMP = Path(__file__).parent
P='http://schemas.openxmlformats.org/presentationml/2006/main'
A='http://schemas.openxmlformats.org/drawingml/2006/main'
R='http://schemas.openxmlformats.org/officeDocument/2006/relationships'
for pre,uri in [('p',P),('a',A),('r',R)]: E.register_namespace(pre,uri)
def tag(ns,n): return '{'+ns+'}'+n
def sub(p,ns,n,**attrs): return E.SubElement(p,tag(ns,n),{k:str(v) for k,v in attrs.items()})
GREEN='086650'; INK='142923'; MUTED='50665D'; LINE='CFDBD3'; PALE='EDF4EF'; GOLD='AC732B'
slides=[]; previews=[]
class Slide:
    def __init__(self,title,subtitle,num,source):
        self.root=E.Element(tag(P,'sld')); common=sub(self.root,P,'cSld'); self.tree=sub(common,P,'spTree')
        nv=sub(self.tree,P,'nvGrpSpPr'); sub(nv,P,'cNvPr',id=1,name=''); sub(nv,P,'cNvGrpSpPr'); sub(nv,P,'nvPr')
        gp=sub(self.tree,P,'grpSpPr'); xf=sub(gp,A,'xfrm')
        for name,attrs in [('off',dict(x=0,y=0)),('ext',dict(cx=0,cy=0)),('chOff',dict(x=0,y=0)),('chExt',dict(cx=0,cy=0))]: sub(xf,A,name,**attrs)
        cm=sub(self.root,P,'clrMapOvr'); sub(cm,A,'masterClrMapping')
        self.im=Image.new('RGB',(1280,720),'white'); self.d=ImageDraw.Draw(self.im)
        self.idx=2
        self.text(70,47,1140,62,title,31.5,INK,True)
        self.text(72,115,1136,37,subtitle,16.5,MUTED)
        self.rect(70,660,1140,1,LINE)
        self.text(72,669,1080,26,'근거  /  '+source,10.5,MUTED)
        self.text(1161,666,48,28,f'{num:02}',12,GREEN,True,'r')
    def shape(self,x,y,w,h,name,fill=None):
        assert x>=0 and y>=0 and x+w<=1281 and y+h<=721,(name,x,y,w,h)
        s=sub(self.tree,P,'sp'); nv=sub(s,P,'nvSpPr'); sub(nv,P,'cNvPr',id=self.idx,name=name); self.idx+=1
        sub(nv,P,'cNvSpPr'); sub(nv,P,'nvPr'); sp=sub(s,P,'spPr'); xf=sub(sp,A,'xfrm')
        sub(xf,A,'off',x=round(x*9525),y=round(y*9525)); sub(xf,A,'ext',cx=round(w*9525),cy=round(h*9525))
        geo=sub(sp,A,'prstGeom',prst='rect'); sub(geo,A,'avLst')
        if fill: sub(sub(sp,A,'solidFill'),A,'srgbClr',val=fill)
        else: sub(sp,A,'noFill')
        sub(sub(sp,A,'ln',w=0),A,'noFill')
        return s
    def rect(self,x,y,w,h,fill):
        self.shape(x,y,w,h,'도형',fill); self.d.rectangle((x,y,x+w,y+h),fill='#'+fill)
    def text(self,x,y,w,h,text,size=18,color=INK,bold=False,align='l'):
        s=self.shape(x,y,w,h,text[:24]); tx=sub(s,P,'txBody'); sub(tx,A,'bodyPr',wrap='square',lIns=0,rIns=0,tIns=0,bIns=0,anchor='ctr'); sub(tx,A,'lstStyle')
        lines=text.split('\n'); fh=size*1.333
        assert len(lines)*fh*1.13 <= h+7,('height',text)
        font=ImageFont.truetype('C:/Windows/Fonts/malgunbd.ttf' if bold else 'C:/Windows/Fonts/malgun.ttf',round(fh))
        for i,line in enumerate(lines):
            p=sub(tx,A,'p'); pp=sub(p,A,'pPr',algn=align); sub(pp,A,'lnSpc').append(E.Element(tag(A,'spcPct'),{'val':'110000'}))
            run=sub(p,A,'r'); rp=sub(run,A,'rPr',lang='ko-KR',sz=round(size*100),b=int(bold)); sub(sub(rp,A,'solidFill'),A,'srgbClr',val=color)
            for n in ['latin','ea','cs']: sub(rp,A,n,typeface='Noto Sans CJK KR')
            sub(run,A,'t').text=line
            length=self.d.textlength(line,font=font)
            # Preview uses a locally available Korean fallback font.
            assert length<=w+12,('width',line,length,w)
            px=x if align=='l' else x+(w-length)/2 if align=='ctr' else x+w-length
            py=y+(h-len(lines)*fh*1.1)/2+i*fh*1.1
            self.d.text((px,py),line,font=font,fill='#'+color)
    def line(self,x1,y1,x2,y2,color=GREEN,width=2):
        # Editable line as a thin rotated rectangle.
        length=math.hypot(x2-x1,y2-y1); s=self.shape(0,0,length,width,'선',color)
        xf=s.find('.//'+tag(A,'xfrm')); xf.set('rot',str(round(math.degrees(math.atan2(y2-y1,x2-x1))*60000)))
        off=xf.find(tag(A,'off')); off.set('x',str(round(((x1+x2)/2-length/2)*9525))); off.set('y',str(round(((y1+y2)/2-width/2)*9525)))
        self.d.line((x1,y1,x2,y2),fill='#'+color,width=width)
    def done(self):
        slides.append(E.tostring(self.root,encoding='utf-8',xml_declaration=True)); previews.append(self.im)

s=Slide('Landing Area와의 수직이격 정의','평면상 중첩 구간에서 항로와 경사 표면 사이의 최소 절대 고도 차를 평가',5,'optimization_path.py · landing_metrics()')
s.text(72,181,570,48,'g = min |z_route - z_area|',26,GREEN,True)
s.text(72,244,558,64,'같은 XY 위치에서 두 고도를 비교하고\n겹치는 전체 구간의 최솟값을 사용합니다.',18)
s.rect(72,337,540,1,LINE)
s.text(72,357,540,35,'평가 범위',18,GREEN,True)
s.text(72,403,558,102,'• 양의 길이로 겹치는 Landing Area마다 계산\n• 자기 FAF의 Landing Area도 포함\n• 미중첩·경계의 점 접촉은 평가에서 제외',16.5)
s.text(72,550,550,67,'위·아래는 절댓값으로 동일하게 평가하며,\n항로가 표면을 가로지르면 g = 0m입니다.',17,MUTED)
s.rect(672,185,536,428,PALE)
s.text(700,201,480,36,'중첩 구간의 측면 개념도',17,GREEN,True)
s.line(720,482,1160,392,MUTED,4); s.text(902,484,280,37,'Landing Area 경사 표면',13,MUTED)
s.line(720,333,1160,357,GREEN,4); s.text(730,277,300,38,'항로 고도 z_route',15,GREEN,True)
s.line(1065,353,1065,411,GOLD,3); s.line(1057,353,1073,353,GOLD,2); s.line(1057,411,1073,411,GOLD,2)
s.text(1080,366,102,38,'고도 차',13,GOLD,True)
s.text(700,549,480,45,'개념도 · 표면의 법선거리가 아닌 수직 고도 차',12.5,MUTED)
s.done()

s=Slide('수직 근접 페널티의 감소 기준','고도 차가 작을수록 비용 증가 · 영역별 비용을 합산하여 후보 점수에 반영',6,'optimization_path.py · Config / clearance_cost / score')
s.rect(72,174,1136,88,PALE)
s.text(94,184,1092,66,'P_vertical = 50,000 × Σ exp(-g / 100)',28,GREEN,True,'ctr')
s.text(72,283,610,35,'100m는 지수함수의 감쇠 척도',19,GREEN,True)
# Editable curve and axes.
x0,y0,xw,yh=117,570,510,220
for val in [0,25000,50000]:
    y=y0-yh*val/50000; s.line(x0,y,x0+xw,y,LINE,1)
    s.text(72,y-12,40,24,f'{val//1000}k',10.5,MUTED,False,'r')
for g in [0,100,200,300]:
    x=x0+xw*g/300; s.text(x-25,y0+8,50,27,str(g),11,MUTED,False,'ctr')
for g in range(0,300,3):
    s.line(x0+xw*g/300,y0-yh*math.exp(-g/100),x0+xw*(g+3)/300,y0-yh*math.exp(-(g+3)/100),GREEN,3)
s.text(265,611,320,29,'최소 절대 고도 차 g (m)',12,MUTED)
s.text(73,319,180,28,'영역당 페널티',11,MUTED)
xs=[727,947]; widths=[220,261]
for j,(gap,pen) in enumerate([('최소 고도 차','영역당 페널티'),('0m','50,000'),('50m','30,327'),('100m','18,394'),('200m','6,767'),('300m','2,489')]):
    y=289+j*44
    s.rect(727,y,481,43,GREEN if j==0 else PALE if j%2 else 'FFFFFF')
    for k,t in enumerate([gap,pen]): s.text(xs[k]+12,y+2,widths[k]-24,39,t,16,'FFFFFF' if j==0 else INK,j==0,'ctr')
s.text(727,573,481,62,'100m에서도 비용은 약 36.8% 남습니다.\n특정 거리에서 0이 되는 구조는 아닙니다.',15,MUTED)
s.done()

s=Slide('후보 점수에 반영되는 방식','수직 근접 비용은 관통 개수·중첩 길이 비용과 함께 더해지는 soft penalty',7,'optimization_path.py · generate() / score · 예시는 설명용 가상 후보')
s.rect(72,177,1136,75,PALE)
s.text(93,187,1094,54,'score = 길이·뱅크 비용 + 관통 비용 + 중첩 비용 + 수직 근접 비용',20,GREEN,True,'ctr')
s.text(72,271,1136,36,'예시 조건  /  Landing Area 1개와 20m 중첩 · 최소 고도 차만 변경',18,INK,True)
rows=[['평가 항목','후보 A · g = 100m','후보 B · g = 200m'],['관통 개수 × 10,000','10,000','10,000'],['중첩 길이 × 10','200','200'],['50,000 × exp(-g / 100)','18,394','6,767'],['관련 비용 합계','28,594','16,967']]
for j,row in enumerate(rows):
    y=327+j*49; s.rect(72,y,1136,48,GREEN if j==0 else PALE if j in [2,4] else 'FFFFFF')
    for k,(x,w) in enumerate([(84,460),(554,314),(878,314)]): s.text(x,y+3,w,42,row[k],16,'FFFFFF' if j==0 else GREEN if j==4 else INK,j in [0,4], 'l' if k==0 else 'ctr')
s.text(72,590,1136,49,'같은 길이·뱅크 비용이라면, 고도 차가 큰 후보 B의 점수가 약 11,627 낮아집니다.',17,GREEN,True)
s.done()

s=Slide('설정값과 적용 단계의 해석','후보 비용의 수직 근접 항과 항로 간 교점 Z 이격 최적화는 서로 다른 평가',8,'optimization_path.py · Config / landing_metrics() · optimization_new_options.py')
s.text(72,179,540,40,'설정값은 Config에 정의',20,GREEN,True)
s.rect(72,236,540,95,PALE)
s.text(90,246,504,74,'clearance_weight = 50,000\nclearance_scale_m = 100m',20,GREEN,True)
s.text(72,352,540,117,'가중치 증가 → 수직 근접 비용의 영향 증가\n감쇠 척도 증가 → 더 먼 거리까지 비용 유지\n여러 영역과 중첩 → 영역별 비용을 합산',16.5)
s.text(72,505,540,108,'현재 코드에 설정된 목적함수 계수입니다.\n100m를 만족하면 통과하는 하드 필터나\n인증된 안전거리라는 의미는 없습니다.',16.5,MUTED)
s.rect(665,181,1,434,LINE)
s.text(710,179,498,40,'선정 과정에서 구분할 점',20,GREEN,True)
s.text(710,245,498,34,'01  후보 필터와 점수는 별도',17,GREEN,True)
s.text(710,287,498,71,'Landing Area 무관통 후보를 우선하는\n필터는 별도로 적용됩니다.',16.5)
s.text(710,389,498,34,'02  측정 대상이 다름',17,GREEN,True)
s.text(710,431,498,71,'이 페널티: 항로 ↔ Landing Area 표면\n교점 Z 이격: 항로 ↔ 다른 항로',16.5)
s.text(710,543,498,71,'기존 자료의 최소 교점 Z 차 13.52m는\n이 페널티의 기준값이 아닙니다.',16.5,MUTED)
s.done()

s=Slide('Hard constraint와 Soft constraint','1차 개별 후보 평가에서 제약 조건과 가중 점수의 역할을 구분',9,'수정본_20261103/코드/optimization_path.py · Config (59행) / generate()')
s.rect(72,180,548,83,PALE)
s.text(94,193,504,56,'Hard constraint · 후보 제외 조건',23,GREEN,True)
s.rect(660,180,548,83,PALE)
s.text(682,193,504,56,'Soft constraint · 가중 점수 평가',23,GREEN,True)
s.text(72,286,548,120,'위반 후보는 평가 대상에서 제외합니다.\n가중치를 곱해 위반을 상쇄하지 않습니다.\n예: 최대 뱅크각·허용 구배 조건',18)
s.text(660,286,548,120,'각 평가 항목에 가중치를 곱해 합산합니다.\n점수가 낮은 후보를 선호합니다.\n예: 관통·중첩·수직 근접·길이·뱅크 비용',18)
s.rect(72,441,1136,1,LINE)
s.text(72,462,1136,39,'가중치 정의 위치',20,GREEN,True)
s.text(72,511,1136,40,'수정본_20261103/코드/optimization_path.py  →  Config 클래스 (59행부터)',18,INK,True)
s.text(72,568,1136,66,'별도 설정 파일이 아닌 코드의 기본값입니다. 실제 실행에 사용된 값은 결과물 summary.json에도 기록됩니다.\nLanding Area 무관통 후보를 우선하는 필터는 가중합 점수와 별도로 적용됩니다.',15.5,MUTED)
s.done()

s=Slide('1차 개별 후보 평가의 가중치','Config 기본값 · 항목별 단위와 정규화 방식에 따라 비용을 계산한 뒤 합산',10,'수정본_20261103/코드/optimization_path.py · Config (76~81행) / score (337~340행)')
rows=[['변수','값','비용에 반영하는 항목'],['crossing_weight','10,000','Landing Area 관통 개수'],['overlap_weight','10','Landing Area 중첩 길이 1m당'],['clearance_weight','50,000','Landing Area와의 수직 근접 페널티'],['length_weight','1','항로 길이 1km당'],['bank_weight','1','최대 허용 뱅크각 대비 비율']]
for j,row in enumerate(rows):
    y=179+j*55; s.rect(72,y,1136,54,GREEN if j==0 else PALE if j%2 else 'FFFFFF')
    for k,(x,w) in enumerate([(88,355),(459,152),(630,560)]):
        s.text(x,y+5,w,44,row[k],17,'FFFFFF' if j==0 else INK,j==0,'ctr' if k==1 else 'l')
s.rect(72,534,1136,101,PALE)
s.text(91,543,1098,37,'수직 근접 항 = clearance_weight × Σ exp(-g / clearance_scale_m)',19,GREEN,True)
s.text(91,586,1098,36,'g: 영역별 최소 절대 고도 차(m)  /  clearance_scale_m = 100m  /  항로 간 교차 가중치는 별도',14.5,MUTED)
s.done()

with zipfile.ZipFile(SRC) as zin:
    entries={n:zin.read(n) for n in zin.namelist()}
    # Retain template masters, theme, dimensions and layout.
    for i,data in enumerate(slides,1):
        entries[f'ppt/slides/slide{i}.xml']=data
        entries[f'ppt/slides/_rels/slide{i}.xml.rels']=b'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideLayout" Target="../slideLayouts/slideLayout6.xml"/></Relationships>'''
    # Remove unreferenced original speaker notes and images to prevent stale content.
    for n in list(entries):
        if n.startswith('ppt/notesSlides/') or n.startswith('ppt/media/'): del entries[n]
    ct=E.fromstring(entries['[Content_Types].xml'])
    for n in list(ct):
        if n.attrib.get('PartName','').startswith('/ppt/notesSlides/'): ct.remove(n)
    relns='http://schemas.openxmlformats.org/package/2006/relationships'
    ctns='http://schemas.openxmlformats.org/package/2006/content-types'
    presentation=E.fromstring(entries['ppt/presentation.xml'])
    ids=presentation.find(tag(P,'sldIdLst'))
    rels=E.fromstring(entries['ppt/_rels/presentation.xml.rels'])
    for i in range(5,len(slides)+1):
        rid=f'rId{20+i}'
        E.SubElement(ids,tag(P,'sldId'),{'id':str(300+i),tag(R,'id'):rid})
        E.SubElement(rels,tag(relns,'Relationship'),{'Id':rid,'Type':R+'/slide','Target':f'slides/slide{i}.xml'})
        E.SubElement(ct,tag(ctns,'Override'),{'PartName':f'/ppt/slides/slide{i}.xml','ContentType':'application/vnd.openxmlformats-officedocument.presentationml.slide+xml'})
    entries['ppt/presentation.xml']=E.tostring(presentation,encoding='utf-8',xml_declaration=True)
    entries['ppt/_rels/presentation.xml.rels']=E.tostring(rels,encoding='utf-8',xml_declaration=True)
    entries['[Content_Types].xml']=E.tostring(ct,encoding='utf-8',xml_declaration=True)
    app=E.fromstring(entries['docProps/app.xml'])
    for el in app:
        if el.tag.endswith('}Slides'): el.text=str(len(slides))
    entries['docProps/app.xml']=E.tostring(app,encoding='utf-8',xml_declaration=True)
    with zipfile.ZipFile(OUT,'w',zipfile.ZIP_DEFLATED) as zout:
        for n,data in entries.items(): zout.writestr(n,data)
for i,im in enumerate(previews,5): im.save(TMP/f'slide_{i:02}.png')
contact=Image.new('RGB',(1280,360*math.ceil(len(previews)/2)),'#dddddd')
for i,im in enumerate(previews): contact.paste(im.resize((640,360)),((i%2)*640,(i//2)*360))
contact.save(TMP/'preview.png')
with zipfile.ZipFile(OUT) as z:
    assert z.testzip() is None
    for n in z.namelist():
        if n.endswith('.xml') or n.endswith('.rels'): E.fromstring(z.read(n))
print(OUT)
print(f'Validated: {len(slides)} slides; editable shapes/text; XML, ZIP and text bounds checked.')
