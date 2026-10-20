# -*- coding: utf-8 -*-
"""用 Tectonic (LaTeX) 渲染全部显示公式为透明背景 PNG"""
import os, subprocess, pymupdf

os.makedirs('/mnt/agents/output/eq_tex', exist_ok=True)
OUT = '/mnt/agents/output/eq/'

EQUATIONS = {
'eq_model':   r'\boxed{\,I_O\ddot{\theta}=N\Bigl(L\sin\theta+\mu h\,\tanh\frac{v_{rel}}{\epsilon}\Bigr)-\frac{mgL}{2}\sin\theta,\qquad I_O=\frac{mL^{2}}{3}\,}',
'eq_constit': r'N=k_c\,(L\cos\theta-h)^{+},\qquad v_{rel}=v-L\cos\theta\,\dot{\theta}',
'eq_tstar':   r'\theta_*=\arccos\frac{h}{L},\qquad N(t_*)=0',
'eq_ic':      r'\theta_0=\arccos\frac{h+N_0/k_c}{L},\qquad \dot{\theta}(0)=\omega_0',
'eq_start':   r'N_0\bigl(L\sin\theta_0+\mu h\bigr)>\frac{mgL}{2}\sin\theta_0\ \ \Longleftrightarrow\ \ N_0\gtrsim\frac{mg\sin\theta_0}{2\bigl(\sin\theta_0+\mu h/L\bigr)}',
'eq_obs':     r'I_{2nd}=\mu\int_{0}^{t_*}N\,|v_{rel}|\,dt,\qquad \ell=\int_{0}^{t_*}|v_{rel}|\,dt,\qquad \mu\langle N\rangle=\frac{I_{2nd}}{\ell}',
'eq_flip':    r'\dot{\theta}_{\mathrm{flip}}=\frac{v}{L\cos\theta}',
'eq_linear':  r'N\approx k_c L\sin\theta_*\,(\theta_*-\theta)',
'eq_energy':  r't_*=\int_{\theta_0}^{\theta_*}\frac{d\theta}{\dot{\theta}(\theta)},\qquad \frac{1}{2}I_O\dot{\theta}^{\,2}=\frac{1}{2}I_O\omega_0^{\,2}+V(\theta_0)-V(\theta)',
'eq_torque':  r'\tau_z=\frac{mgL}{2}\sin\theta-NL\sin\theta-\mu h N\,\tanh\frac{v_{rel}}{\epsilon},\qquad \ddot{\theta}=-\tau_z/I_O',
'eq_deposit': r'q(t)=\alpha\mu N(t)\,|v_{rel}(t)|',
'eq_ximap':   r'\xi(t)=\int_{0}^{t}v_{rel}(\tau)\,d\tau,\qquad \ell=\int_{0}^{t_*}|v_{rel}|\,dt',
'eq_foldsum': r'I(\xi)=\sum_{\{i\,:\,\xi(t_i)=\xi\}}\alpha\mu\,N(t_i)',
}

TEMPLATE = r'''\documentclass[border=3pt]{standalone}
\usepackage{amsmath,amssymb}
\begin{document}
$\displaystyle %s$
\end{document}
'''

for name, tex in EQUATIONS.items():
    d = f'/mnt/agents/output/eq_tex/{name}'
    os.makedirs(d, exist_ok=True)
    with open(f'{d}/{name}.tex', 'w') as f:
        f.write(TEMPLATE % tex)
    subprocess.run(['tectonic', f'{name}.tex'], cwd=d,
                   capture_output=True, timeout=240)
    doc = pymupdf.open(f'{d}/{name}.pdf')
    pix = doc[0].get_pixmap(dpi=300, alpha=True)
    pix.save(OUT + name + '.png')
    print(name, f'{pix.width}x{pix.height}')
print('all done')
