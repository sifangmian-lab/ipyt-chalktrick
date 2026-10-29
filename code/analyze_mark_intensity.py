"""Explore raw-image chalk brightness profiles; no force/absolute deposit inversion.

python code/analyze_mark_intensity.py --photo1 IMAGE --photo2 IMAGE --out DIRECTORY
Uses the existing manual centerline and candidate segments from paired-audio analysis.
"""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from PIL import Image
from scipy import ndimage, signal
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties


def write_csv(p, header, rows):
    with p.open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.writer(f); w.writerow(header); w.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for s in ['photo1', 'photo2', 'out']:
        parser.add_argument('--'+s, type=Path, required=True)
    a = parser.parse_args(); a.out.mkdir(parents=True, exist_ok=True)
    repo = Path(__file__).resolve().parents[1]
    prior = repo/'experiments/E5E6_20261029_audio_marks'
    font = Path('C:/Windows/Fonts/msyh.ttc')
    if font.exists():
        plt.rcParams['font.family'] = FontProperties(fname=str(font)).get_name()
    plt.rcParams['axes.unicode_minus'] = False
    cfg = {
        1: ([150,300,600,900,1200,1500,1800],[693,678,646,619,614,605,614],4,[12,25,36,44]),
        2: ([0,300,600,900,1200,1500,1700],[606,595,604,611,609,604,607],3,[8,11,15,20])}
    summary = {'status': 'exploratory optical profiles; no radiometric or geometric calibration',
               'direction': 'left_to_right', 'coordinates': 'original image pixels',
               'photos': {}}
    selected = []
    for k in [1,2]:
        im = np.asarray(Image.open(getattr(a,'photo'+str(k))).convert('RGB'),dtype=float)
        gray = im@np.array([.2126,.7152,.0722])
        scale = im.shape[1]/1824
        xg,yg,deg,selection = cfg[k]
        x = np.arange(gray.shape[1]); y = np.polyval(np.polyfit(xg,yg,deg),x/scale)*scale
        yy = np.arange(-round(65*scale),round(65*scale)+1)
        coords = np.array([y[None,:]+yy[:,None],np.broadcast_to(x,(len(yy),len(x)))])
        band = ndimage.map_coordinates(gray,coords,order=1)
        upper = (yy<=-45*scale)&(yy>=-60*scale)
        lower = (yy>=45*scale)&(yy<=60*scale)
        bu = ndimage.gaussian_filter1d(np.median(band[upper],axis=0),3)
        bl = ndimage.gaussian_filter1d(np.median(band[lower],axis=0),3)
        bg = bu[None,:]+(bl-bu)[None,:]*(yy[:,None]-yy[upper].mean())/(yy[lower].mean()-yy[upper].mean())
        contrast = band-bg
        strip = abs(yy)<=38*scale
        c = contrast[strip]
        J = ndimage.gaussian_filter1d(c.sum(axis=0),1)
        widths = {v:ndimage.gaussian_filter1d((c>v).sum(axis=0).astype(float),1) for v in [15,25,35]}
        W = widths[25]
        # Mean brightness among mask-positive pixels, separate from integrated brightness.
        density = ndimage.gaussian_filter1d(np.where(c>25,c,0).sum(axis=0),1)/np.maximum(W,1)
        with (prior/f'photo_{k}_candidate_geometry.csv').open(encoding='utf-8') as f:
            marks = list(csv.DictReader(f))
        allmetrics=[]
        for m in marks:
            idx=int(m['candidate_index'])
            if m['status']!='interior_geometry':continue
            center=float(m['center_x_px'])*scale
            left=int(float(m['left_px'])*scale);right=int(float(m['right_px'])*scale)
            pad=round(10*scale)
            lo=max(0,left-pad);hi=min(len(x),right+pad+1)
            xx=x[lo:hi];jj=J[lo:hi].copy()
            floor=float(np.median(np.r_[jj[:max(2,pad//2)],jj[-max(2,pad//2):]]))
            jj=np.maximum(jj-floor,0)
            peak=jj.max(); norm=jj/peak
            # Use the main connected component; isolated dust is not a trailing edge.
            mask=ndimage.binary_closing(norm>.1,structure=np.ones(3))
            labs,_=ndimage.label(mask)
            component=labs[np.argmax(norm)]
            ids=np.flatnonzero(labs==component) if component else np.array([])
            if len(ids)<3:continue
            L,R=int(ids[0]),int(ids[-1]); length=R-L+1
            # Spatial 10-90 edges are descriptive and include contact footprint and optical blur.
            ids90=np.flatnonzero(norm>=.9)
            rise=float(ids90[0]-L);fall=float(R-ids90[-1])
            support=(np.arange(len(jj))>=L)&(np.arange(len(jj))<=R)
            occupied=(W[lo:hi]>.2*W[lo:hi].max())&support
            corr=float(np.corrcoef(jj[occupied],W[lo:hi][occupied])[0,1]) if occupied.sum()>3 else None
            density_cv=float(np.std(density[lo:hi][occupied])/np.mean(density[lo:hi][occupied]))
            brightpix=(band[strip,lo:hi][c[:,lo:hi]>25])
            sat=float(np.mean(brightpix>=250)) if len(brightpix) else 0
            saturated=float(np.mean(brightpix>=254.5)) if len(brightpix) else 0
            correlations={}
            for threshold,ww in widths.items():
                active=(ww[lo:hi]>.2*ww[lo:hi].max())&support
                correlations[str(threshold)]=float(np.corrcoef(jj[active],ww[lo:hi][active])[0,1]) if active.sum()>3 else None
            area=float(jj.sum())
            peak_u=float((np.argmax(jj)-L)/max(R-L,1))
            met=dict(mark=idx,x_center_px=center,length_10pct_px=length,
                     local_profile_floor=floor,
                     rise_10_90_px=rise,fall_90_10_px=fall,
                     fall_fraction=fall/length,peak_fraction=peak_u,
                     width_intensity_correlation=corr,
                     masked_brightness_cv=density_cv,bright_pixel_near_saturation_fraction=sat,
                     bright_pixel_at_ceiling_fraction=saturated,
                     width_intensity_correlations_by_threshold=correlations,
                     relative_integrated_brightness=area)
            allmetrics.append(met)
            if idx in selection:
                selected.append(dict(photo=k,mark=idx,x=xx-center,
                                     patch=band[strip,lo:hi],J=norm,
                                     W=W[lo:hi]/max(W[lo:hi].max(),1),
                                     density=np.where(occupied,density[lo:hi],np.nan),L=L,R=R,metrics=met))
                write_csv(a.out/f'photo{k}_mark{idx}_profile.csv',
                          ['x_original_px','raw_background_subtracted_integrated_brightness',
                           'local_baseline_adjusted_nonnegative_brightness',
                           'normalized_integrated_brightness','width_threshold25_px',
                                     'mean_contrast_above_threshold25'],
                          zip(xx,J[lo:hi],jj,norm,W[lo:hi],density[lo:hi]))
        summary['photos'][str(k)]={'metrics':allmetrics,'selected_marks':selection,
                                    'background': 'linear across y, side strips 45-60 resized pixels from centerline',
                                    'integration_halfwidth_resized_px':38,
                                    'original_to_1824_scale':scale}
        write_csv(a.out/f'photo{k}_all_metrics.csv',list(allmetrics[0]),
                  [[m[f] for f in allmetrics[0]] for m in allmetrics])

    fig,axs=plt.subplots(2,4,figsize=(14,7),layout='constrained',gridspec_kw={'height_ratios':[1.15,1]})
    chosen=[s for s in selected if s['photo']==2]
    for col,s in enumerate(chosen):
        axs[0,col].imshow(s['patch'],cmap='gray',vmin=0,vmax=255,
                          extent=(s['x'][0],s['x'][-1],-38,38),aspect='auto')
        axs[0,col].set(title=f'图2，第{s["mark"]}个候选笔痕',xlabel='x（原图像素，书写方向 →）')
        axs[0,col].set_yticks([])
        axs[1,col].plot(s['x'],s['J'],label='按 y 汇总的亮度',color='#175d91')
        axs[1,col].plot(s['x'],s['W'],'--',label='亮区有效宽度',color='#cf6b27')
        axs[1,col].set(xlabel='x（原图像素）',ylim=(-.03,1.12),ylabel='各自归一化')
        axs[1,col].legend(fontsize=8)
    fig.suptitle('沿 x 的亮度峰与笔痕形状同时变化；光学亮度尚未标定为沉积量',fontsize=13)
    fig.savefig(a.out/'intensity_width_comparison.png',dpi=150);plt.close(fig)
    fig,axs=plt.subplots(2,4,figsize=(14,7),layout='constrained')
    for row,k in enumerate([1,2]):
        for col,s in enumerate([z for z in selected if z['photo']==k]):
            axs[row,col].plot(s['x'],s['J'],label='积分亮度/峰值')
            axs[row,col].plot(s['x'],s['density']/max(np.nanmax(s['density']),1),label='亮区平均对比度/峰值',alpha=.7)
            axs[row,col].set(title=f'图{k} 笔痕{s["mark"]}',xlabel='x（原图像素）',ylim=(0,1.1))
            axs[row,col].legend(fontsize=7)
    fig.savefig(a.out/'selected_profiles.png',dpi=150);plt.close(fig)
    (a.out/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps([s['metrics']|{'photo':s['photo']} for s in selected],indent=2,ensure_ascii=False))


if __name__=='__main__':main()
