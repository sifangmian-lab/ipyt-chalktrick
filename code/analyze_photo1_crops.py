"""Crop six photo-1 marks without rescaling pixels; compare fixed-height profiles."""
import argparse,csv,json
from pathlib import Path
import numpy as np
from PIL import Image
from scipy import ndimage
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

def write_csv(path,header,rows):
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.writer(f);w.writerow(header);w.writerows(rows)

def width_at_half(p):
    mask=p>=.5
    labs,_=ndimage.label(mask)
    lab=labs[np.argmax(p)]
    ids=np.flatnonzero(labs==lab) if lab else []
    return int(len(ids))

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--photo',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    repo=Path(__file__).resolve().parents[1]
    with (repo/'experiments/E5E6_20261029_audio_marks/photo_1_candidate_geometry.csv').open(encoding='utf-8') as f:
        marks={int(m['candidate_index']):m for m in csv.DictReader(f)}
    im=Image.open(a.photo).convert('RGB');scale=im.width/1824
    centerline=np.polyfit([150,300,600,900,1200,1500,1800],[693,678,646,619,614,605,614],4)
    selected=[12,21,30,36,44,51];all_data=[]
    for idx in selected:
        center=int(round(float(marks[idx]['center_x_px'])*scale))
        yc=int(round(np.polyval(centerline,center/scale)*scale))
        box=(center-48,yc-80,center+48,yc+80)
        patch=im.crop(box);patch.save(a.out/f'photo1_mark{idx}_raw.png')
        rgb=np.asarray(patch,dtype=float);G=rgb@np.array([.2126,.7152,.0722])
        xx=np.arange(G.shape[1])-48;yy=np.arange(G.shape[0])-80
        top=ndimage.gaussian_filter1d(np.median(G[:20],axis=0),3)
        bot=ndimage.gaussian_filter1d(np.median(G[-20:],axis=0),3)
        B=top[None,:]+(bot-top)[None,:]*(np.arange(len(yy))[:,None]-9.5)/140
        D=G-B
        # Select main bright connected footprint near the previously located mark.
        labs,n=ndimage.label(ndimage.gaussian_filter(D,.6)>25)
        scores=[]
        for lab in range(1,n+1):
            rr,cc=np.where(labs==lab)
            if len(rr)<5:continue
            score=np.maximum(D[rr,cc],0).sum()/(1+(np.mean(cc)-48)**2/100+(np.mean(rr)-80)**2/1600)
            scores.append((score,lab))
        lab=max(scores)[1];rr,cc=np.where(labs==lab)
        y0,y1=int(rr.min()),int(rr.max())
        # Integrate a fixed window enclosing the main footprint; keep signed values.
        lo=max(20,y0-4);hi=min(140,y1+5)
        J=ndimage.gaussian_filter1d(D[lo:hi].sum(axis=0),.6)
        floor=float(np.median(np.r_[J[:15],J[-15:]]));J=J-floor
        normJ=J/J.max()
        band_y=[int(round(y0+f*(y1-y0))) for f in [.25,.5,.75]]
        profiles=[];band_stats=[]
        for row in band_y:
            p=ndimage.gaussian_filter1d(D[row-2:row+3].mean(axis=0),.6)
            p-=np.median(np.r_[p[:15],p[-15:]])
            pn=p/p.max()
            profiles.append(pn)
            band_stats.append({'y_original_px':box[1]+row,'peak_x_relative_px':int(xx[np.argmax(pn)]),
                               'fwhm_px':width_at_half(pn),'contrast_peak_gray':float(p.max())})
        sensitivity={}
        for band_height in [3,5,9]:
            peaks=[]
            for row in band_y:
                half=band_height//2
                p=ndimage.gaussian_filter1d(D[row-half:row+half+1].mean(axis=0),.6)
                peaks.append(int(xx[np.argmax(p)]))
            sensitivity[str(band_height)]={'peak_x_relative_px':peaks,'spread_px':int(np.ptp(peaks))}
        weight=np.maximum(D[rr,cc],0)
        mx=np.average(cc,weights=weight);my=np.average(rr,weights=weight)
        slope=float(np.sum(weight*(rr-my)*(cc-mx))/np.sum(weight*(rr-my)**2))
        sat=float(np.mean(G[rr,cc]>=254.5))
        met={'mark':idx,'original_crop_box':box,'main_footprint_bbox_relative':[int(cc.min()),y0,int(cc.max()),y1],
             'integrated_fwhm_px':width_at_half(normJ),'integrated_peak_x_relative_px':int(xx[np.argmax(normJ)]),
             'top_middle_bottom':band_stats,
             'band_peak_x_spread_px':int(np.ptp([q['peak_x_relative_px'] for q in band_stats])),
             'band_height_sensitivity':sensitivity,
             'brightness_weighted_dx_dy':slope,'main_component_at_gray_ceiling_fraction':sat,
             'integrated_profile_local_floor':floor}
        write_csv(a.out/f'photo1_mark{idx}_profiles.csv',
                  ['x_original_px','signed_integrated_contrast','integrated_normalized','top_5rows_normalized','middle_5rows_normalized','bottom_5rows_normalized'],
                  zip(xx+center,J,normJ,*profiles))
        all_data.append({'idx':idx,'box':box,'patch':patch,'xx':xx,'yy':yy,'J':normJ,
                         'profiles':profiles,'band_y':band_y,'metrics':met})

    font=Path('C:/Windows/Fonts/msyh.ttc')
    if font.exists():plt.rcParams['font.family']=FontProperties(fname=str(font)).get_name()
    plt.rcParams['axes.unicode_minus']=False
    fig,axs=plt.subplots(2,6,figsize=(17,8),layout='constrained',gridspec_kw={'height_ratios':[1.4,1]})
    colors=['#0072b2','#009e73','#d55e00'];labels=['上部5行','中部5行','下部5行']
    for i,d in enumerate(all_data):
        axs[0,i].imshow(d['patch'],interpolation='nearest',extent=(-48,48,80,-80),aspect='equal')
        for row,color in zip(d['band_y'],colors):axs[0,i].axhline(row-80,color=color,lw=.8,alpha=.8)
        axs[0,i].set(title=f'图1 编号{d["idx"]}',xlabel='x（原图像素）',ylabel='y（原图像素）' if i==0 else '')
        axs[0,i].set_yticks([-60,0,60])
        axs[1,i].plot(d['xx'],d['J'],color='black',lw=2,label='整痕按y积分')
        for p,color,label in zip(d['profiles'],colors,labels):axs[1,i].plot(d['xx'],p,color=color,lw=1,label=label)
        axs[1,i].set(xlim=(-30,30),ylim=(-.15,1.15),xlabel='x（原图像素）',ylabel='各自归一化亮度' if i==0 else '')
        axs[1,i].axhline(0,color='gray',lw=.5)
        axs[1,i].set_title(f'三条带峰位跨度 {d["metrics"]["band_peak_x_spread_px"]} px',fontsize=10)
    handles,labels=axs[1,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='outside lower center',ncol=4,frameon=False)
    fig.suptitle('图1：原比例裁片与上／中／下亮度剖面（书写方向 x →）\n裁片保持原像素与长宽比；放大显示不增加分辨率，曲线保留负值',fontsize=14)
    fig.savefig(a.out/'photo1_six_marks.png',dpi=150);plt.close(fig)
    fig,axs=plt.subplots(2,3,figsize=(10,9),layout='constrained')
    for ax,d in zip(axs.ravel(),all_data):
        ax.imshow(d['patch'],interpolation='nearest',aspect='equal')
        ax.set(title=f'图1 编号{d["idx"]}',xticks=[],yticks=[])
    fig.suptitle('图1的六个原始裁片｜书写方向：从左向右\n原图像素直接裁取；统一96×160像素，保持长宽比',fontsize=13)
    fig.savefig(a.out/'photo1_raw_crop_montage.png',dpi=150);plt.close(fig)
    # Overview retains the complete trace so the sample selection is auditable.
    fig,ax=plt.subplots(figsize=(15,3.5),layout='constrained')
    ax.imshow(im);ax.set_ylim(1700,1250);ax.set_xlim(250,4096)
    for d in all_data:
        left,top,right,bottom=d['box']
        ax.add_patch(plt.Rectangle((left,top),right-left,bottom-top,fill=False,color='#e31a1c',lw=1.3))
        ax.text((left+right)/2,top-10,str(d['idx']),color='red',ha='center')
    ax.set(title='图1取样位置：12、21、30、36、44、51（候选编号）',xlabel='原图 x 像素',ylabel='原图 y 像素')
    fig.savefig(a.out/'photo1_sample_locations.png',dpi=150);plt.close(fig)
    result={'source':a.photo.name,'source_dimensions':list(im.size),'selection':selected,
            'method':'unrotated raw RGB crops; same 5-row bands at 25/50/75% main footprint height; signed profiles; Gaussian sigma=0.6 original px',
            'metrics':[d['metrics'] for d in all_data]}
    (a.out/'summary.json').write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps(result['metrics'],indent=2,ensure_ascii=False))

if __name__=='__main__':main()
