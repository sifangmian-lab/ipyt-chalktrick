# -*- coding: utf-8 -*-
"""按需子集化 + OTF->TTF 转换，供 ReportLab 嵌入"""
import os, time
from fontTools.ttLib import TTFont, newTable
from fontTools.pens.cu2quPen import Cu2QuPen
from fontTools.pens.ttGlyphPen import TTGlyphPen
from fontTools.subset import Subsetter, Options

FD = '/usr/share/fonts/opentype/noto/'
OUT = '/mnt/agents/output/fonts/'
os.makedirs(OUT, exist_ok=True)

# 收集正文所需全部字符：手工指定常用汉字 + ASCII + 全角标点
EXTRA = '，。、；：（）—％·《》〈〉「」！？“”‘’ⅠⅡⅢⅣⅤⅥⅦⅧⅨⅩ₀₁₂₃₄₅₆₇₈₉⁻×－～℃'
ASCII = ''.join(chr(c) for c in range(32, 127))

def build_text():
    import re
    text = EXTRA + ASCII
    # 提取两个源文件中的全部非 ASCII 字符（汉字、全角标点、希腊字母、数学符号等）
    src = open('/mnt/agents/output/build_pdf.py', encoding='utf-8').read()
    code = open('/mnt/agents/output/粉笔_N下降段_完整模型.py', encoding='utf-8').read()
    code2 = open('/mnt/agents/output/划痕相对强度_时空依赖.py', encoding='utf-8').read()
    for s in (src, code, code2):
        text += ''.join(re.findall(r'[^\x00-\x7f]', s))
    return ''.join(sorted(set(text)))

TEXT = build_text()
print(f'字符集大小: {len(TEXT)}')

def convert(src_ttc, index, out_name):
    t0 = time.time()
    font = TTFont(src_ttc, fontNumber=index, lazy=True)
    # 子集化
    opt = Options()
    opt.name_IDs = ['*']
    opt.notdef_outline = True
    opt.recalc_bounds = True
    opt.drop_tables += ['vhea', 'vmtx', 'VORG', 'GSUB', 'GPOS', 'GDEF', 'BASE']
    ss = Subsetter(options=opt)
    ss.populate(text=TEXT)
    ss.subset(font)
    # OTF -> TTF
    glyphOrder = font.getGlyphOrder()
    glyphSet = font.getGlyphSet()
    glyf_glyphs = {}
    for name in glyphOrder:
        glyph = glyphSet[name]
        pen = TTGlyphPen(glyphSet)
        cu2qu = Cu2QuPen(pen, 1.0, reverse_direction=True)
        glyph.draw(cu2qu)
        glyf_glyphs[name] = pen.glyph()
    glyf = newTable('glyf')
    glyf.glyphOrder = glyphOrder
    glyf.glyphs = glyf_glyphs
    font['glyf'] = glyf
    font['loca'] = newTable('loca')
    maxp = font['maxp']
    maxp.tableVersion = 0x00010000
    for attr in ('maxPoints maxContours maxCompositePoints maxCompositeContours '
                 'maxTwilightPoints maxStorage maxFunctionDefs '
                 'maxInstructionDefs maxStackElements maxSizeOfInstructions '
                 'maxComponentElements maxComponentDepth').split():
        setattr(maxp, attr, 0)
    maxp.maxZones = 1
    post = font['post']
    post.formatType = 3.0
    post.extraNames = []
    post.mapping = {}
    post.glyphOrder = glyphOrder
    for tag in ('CFF ', 'VORG'):
        if tag in font:
            del font[tag]
    font.sfntVersion = '\x00\x01\x00\x00'
    font.save(OUT + out_name)
    print(f'{out_name}: {os.path.getsize(OUT+out_name)/1e6:.2f} MB, {time.time()-t0:.0f}s')

convert(FD + 'NotoSerifCJK-Regular.ttc', 2, 'SongSC.ttf')
convert(FD + 'NotoSerifCJK-Bold.ttc', 2, 'SongSC-Bold.ttf')
convert(FD + 'NotoSansCJK-Bold.ttc', 2, 'HeiSC-Bold.ttf')
convert(FD + 'NotoSansCJK-Regular.ttc', 7, 'MonoSC.ttf')
print('done')
