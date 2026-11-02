"""S8: forward-map v2 deposition through a finite footprint and compare photo 1.

Image-scale spans and footprint parameters are scenarios, not calibrated values.
Only scalar brightness gain and x translation are fitted to the measured curves.
No dynamical parameter, force history, or source scale is inferred from a photo.
"""
from pathlib import Path
import argparse,csv,json
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from scipy import ndimage
from observation_model import footprint_kernel,profile_to_pixel_mass,render_mass,camera_response,camera_image_on_crop
import v2_model as vm

def write_csv(path,header,rows):
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.writer(f);w.writerow(header);w.writerows(rows)

def observed_crop(path):
    rgb=np.array(Image.open(path).convert('RGB'),dtype=float)
    g=rgb@np.array([.2126,.7152,.0722])
    top=ndimage.gaussian_filter1d(np.median(g[:20],axis=0),3)
    bottom=ndimage.gaussian_filter1d(np.median(g[-20:],axis=0),3)
    bg=top[None,:]+(bottom-top)[None,:]*(np.arange(g.shape[0])[:,None]-9.5)/140
    d=g-bg
    j=d[20:140].sum(axis=0)
    floor=float(np.median(np.r_[j[:15],j[-15:]]))
    return rgb,d,j-floor

def register_projection(x,field,observed):
    """Fit translation and positive gain only; keep signed observation noise."""
    coords=np.arange(len(observed))-len(observed)//2
    window=abs(coords)<=30
    base=field.sum(axis=0)
    best=None
    for shift in np.arange(-15,15.001,.25):
        predicted=np.interp(coords-shift,x,base,left=0,right=0)
        gain=max(0,float(observed[window]@predicted[window]/max(predicted[window]@predicted[window],1e-30)))
        error=float(np.sqrt(np.mean((observed[window]-gain*predicted[window])**2))/max(observed.max(),1e-12))
        if best is None or error<best['nrmse']:
            best=dict(x_shift_px=float(shift),gain=gain,nrmse=error,predicted=gain*predicted)
    return best

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    repo=Path(__file__).resolve().parents[1]
    prior=repo/'experiments/E5E6_20261029_intensity/photo1_examples'
    meta=json.loads((prior/'summary.json').read_text(encoding='utf-8'))
    metrics={m['mark']:m for m in meta['metrics']}
    s4=repo/'simulations/S4_v2正交弹簧基线/S4_基线A_强度剖面.csv'
    data=np.genfromtxt(s4,delimiter=',',skip_header=1,encoding='utf-8-sig')
    sources={'spring-slider':(data[:,0]*1e-3,data[:,1])}
    p=dict(vm.DEFAULT,friction='lugre',sigma1=30.,sigma2=5.)
    info,tr=vm.simulate(p,t_max=.004)
    xc,I,met=vm.dash_profile(tr,nbins=120,sigma_bins=1.)
    sources['viscous']=(xc,I)
    for tag,(x,i) in sources.items():
        write_csv(a.out/f'source_{tag}.csv',['xi_m','I_J_per_m'],zip(x,i))
    # Only slope and height are suggested by mark 44; sigma/PSF are declared scenarios.
    geometry=dict(half_height_px=27,sigma_x_px=1.5,tilt_dx_dy=.31)
    kernels={'point':np.ones((1,1)),
             'upright':footprint_kernel(**(geometry|{'tilt_dx_dy':0})),
             'tilted':footprint_kernel(**geometry)}
    spans=[2.,4.,8.,16.];psf=.7
    observations={m:observed_crop(prior/f'photo1_mark{m}_raw.png') for m in [44,51]}
    rows=[];examples={};mass_errors=[]
    for tag,(xx,ii) in sources.items():
        for span in spans:
            grid,mass=profile_to_pixel_mass(xx,ii,span)
            original_mass=float(mass.sum())
            mass=mass/original_mass  # Shape-only forward experiments; absolute mass stays in source CSV.
            for name,kernel in kernels.items():
                x,y,field=render_mass(mass,kernel,grid[0],psf)
                mass_errors.append(abs(float(field.sum())-1))
                for mark,(_,_,obs) in observations.items():
                    fit=register_projection(x,field,obs)
                    rows.append(dict(source=tag,footprint=name,source_span_px=span,mark=mark,
                                     x_shift_px=fit['x_shift_px'],brightness_gain=fit['gain'],
                                     nrmse=fit['nrmse'],model_total_proxy_before_unit_normalization=original_mass))
                    if tag=='spring-slider' and span==4:
                        examples[(mark,name)]=(x,y,field,fit)
    write_csv(a.out/'comparison_scores.csv',list(rows[0]),[[r[k] for k in rows[0]] for r in rows])
    font=Path('C:/Windows/Fonts/msyh.ttc')
    if font.exists():plt.rcParams['font.family']=FontProperties(fname=str(font)).get_name()
    plt.rcParams['axes.unicode_minus']=False
    labels={'point':'点足迹','upright':'直立有限足迹','tilted':'倾斜有限足迹'}
    colors={'point':'#777777','upright':'#0072b2','tilted':'#d55e00'}
    fig,axs=plt.subplots(2,3,figsize=(14,8),layout='constrained')
    for tag,(xx,ii) in sources.items():
        u=(xx-xx[0])/(xx[-1]-xx[0])
        axs[0,0].plot(u,ii/np.median(ii),label=tag)
    axs[0,0].set(title='输入：现有v2单事件沉积剖面',xlabel='各自归一化空间位置',ylabel='I / 各自中位数')
    axs[0,0].legend(fontsize=9)
    axs[0,1].imshow(kernels['tilted'],origin='lower',aspect='equal',cmap='magma')
    axs[0,1].set(title='有限足迹 K：单位总量，倾斜脊线',xlabel='足迹网格 x',ylabel='足迹网格 y')
    x,y,field,fit=examples[(44,'tilted')]
    crop=camera_image_on_crop(x,y,field,x_shift=fit['x_shift_px'],
                              y_center=np.mean(metrics[44]['main_footprint_bbox_relative'][1::2]))
    axs[0,2].imshow(crop,cmap='gray',origin='upper',aspect='equal')
    axs[0,2].set(title='正向生成二维沉积形状（非原图）',xlabel='裁片像素 x',ylabel='裁片像素 y')
    for col,mark in enumerate([44,51]):
        obs=observations[mark][2];coords=np.arange(len(obs))-len(obs)//2
        axs[1,col].plot(coords,obs/obs.max(),color='black',lw=2,label=f'图1 编号{mark}')
        allcols=[coords,obs]
        for name in kernels:
            xx,yy,dd,ft=examples[(mark,name)]
            axs[1,col].plot(coords,ft['predicted']/obs.max(),label=labels[name],color=colors[name])
            allcols.append(ft['predicted'])
        write_csv(a.out/f'mark{mark}_comparison.csv',
                  ['x_relative_px','observed_signed_contrast','point','upright','tilted'],zip(*allcols))
        axs[1,col].set(title=f'编号{mark}：源跨度固定为4 px',xlim=(-30,30),ylim=(-.15,1.2),
                      xlabel='原图 x 像素（相对裁片中心）',ylabel='以实测峰值归一化')
        axs[1,col].legend(fontsize=8)
    for mark in [44,51]:
        for name in ['point','tilted']:
            rr=[r for r in rows if r['source']=='spring-slider' and r['mark']==mark and r['footprint']==name]
            axs[1,2].plot([r['source_span_px'] for r in rr],[r['nrmse'] for r in rr],
                          'o-' if mark==44 else 's--',color=colors[name],label=f'{mark} {labels[name]}')
    axs[1,2].set(title='未标定空间尺度的情景敏感性',xlabel='指定源轨迹跨度 / px',ylabel='配准后 NRMSE（越小越接近）')
    axs[1,2].legend(fontsize=8)
    fig.suptitle('S8：动力学沉积 → 有限二维足迹 → 成像 → x向剖面\n图像对照只拟合平移与亮度比例；足迹与像素尺度尚未独立标定',fontsize=14)
    fig.savefig(a.out/'forward_comparison.png',dpi=150);plt.close(fig)

    # Show camera nonlinearity order and loss of distinctions between source laws.
    fig,axs=plt.subplots(1,3,figsize=(14,4),layout='constrained')
    grid,mass=profile_to_pixel_mass(*sources['spring-slider'],4)
    mass/=mass.sum()
    x,y,field=render_mass(mass,kernels['tilted'],grid[0],psf)
    gain=1/field.max()
    camera_rows=[]
    for mode,fac in [('linear',1),('saturating',1),('clipped',3)]:
        image=camera_response(field,gain=gain*fac,mode=mode)
        jj=image.sum(axis=0)
        axs[0].plot(x,jj/jj.max(),label=f'{mode}, gain×{fac}')
        camera_rows.append((mode,fac,jj))
    axs[0].set(title='逐像素相机响应后再按y积分',xlabel='x / px',ylabel='各自峰值归一化');axs[0].legend(fontsize=8)
    distinctions=[]
    for col,span in enumerate([2.,16.]):
        results=[]
        for tag,(xx,ii) in sources.items():
            grid,mass=profile_to_pixel_mass(xx,ii,span);mass/=mass.sum()
            xp,yp,dd=render_mass(mass,kernels['tilted'],grid[0],psf)
            jj=dd.sum(axis=0);jj/=jj.max()
            axs[col+1].plot(xp,jj,label=tag)
            results.append((xp,jj))
        ref=np.arange(-40,41,.25)
        aa=np.interp(ref,*results[0],left=0,right=0);bb=np.interp(ref,*results[1],left=0,right=0)
        diff=float(np.linalg.norm(aa-bb)/np.linalg.norm(aa))
        distinctions.append(dict(source_span_px=span,relative_L2_between_peak_normalized_profiles=diff))
        axs[col+1].set(title=f'相同足迹，源跨度{span:g}px\n两种源形状相对L2差异 {diff:.3f}',xlabel='x / px')
        axs[col+1].legend(fontsize=8)
    fig.savefig(a.out/'observation_sensitivity.png',dpi=150);plt.close(fig)
    result=dict(geometry_scenario=geometry,psf_sigma_px=psf,source_span_scenarios_px=spans,
                max_mass_conservation_error=max(mass_errors),
                comparison_fit='positive scalar gain and x translation only, x in [-30,30]',
                calibration_status='pixel scale, footprint and camera uncalibrated; NOT a force inverse fit',
                geometry_provenance='height and tilt suggested by observed mark44 footprint; sigma/PSF chosen',
                comparison_status='44 is exploratory self-comparison; 51 shares fixed geometry but is not independent physical validation',
                viscous_simulation_parameters=p,viscous_source_metrics=met,
                distinctions=distinctions,
                fixed_span4_scores=[r for r in rows if r['source']=='spring-slider' and r['source_span_px']==4])
    (a.out/'summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
