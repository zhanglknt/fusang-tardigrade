import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np, os

fig, ax = plt.subplots(figsize=(10, 4), facecolor='white')
ax.set_xlim(0, 10); ax.set_ylim(0, 4); ax.axis('off')

# Title
ax.text(5, 3.85, 'IMMI: Information-Matched Multi-Level Inference Framework',
         ha='center', va='top', fontsize=14, fontweight='bold')

# Levels
levels = [
    (0.2, 2.8, '#E8F4FD', '#2196F3', 'Level 0', 'k-mer extraction', 'O(n·L) | no alignment'),
    (2.5, 2.8, '#FFF3E0', '#FF9800', 'Level 1', 'Cosine+NJ', 'O(n²) | backbone'),
    (4.8, 2.8, '#E8F5E9', '#4CAF50', 'Level 2', 'RF classifier', 'AUC=0.99 | escalate?'),
    (7.1, 2.8, '#FFEBEE', '#F44336', 'Level 3', 'MSA+ML refine', 'O(m²L²) | subtrees'),
]
for (x, y, fc, ec, title, subtitle, note) in levels:
    ax.add_patch(FancyBboxPatch((x, y), 1.8, 0.6, boxstyle="round,pad=0.05",
                                  facecolor=fc, edgecolor=ec, linewidth=1.5))
    ax.text(x+0.9, y+0.38, title, ha='center', va='center', fontsize=9, fontweight='bold')
    ax.text(x+0.9, y+0.22, subtitle, ha='center', va='center', fontsize=7)
    ax.text(x+0.9, y+0.08, note, ha='center', va='center', fontsize=5.5, color='#555')

# Arrows
for (x1,y1,x2,y2) in [(2.0,3.1,2.5,3.1),(4.3,3.1,4.8,3.1),(6.6,3.1,7.1,3.1)]:
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2), arrowstyle='->',
                                    mutation_scale=15, linewidth=1.5, color='#333',
                                    linestyle='--' if x1>2 else '-'))

# Info gradient
ax.add_patch(mpatches.Rectangle((0.2,2.5),8.6,0.06, facecolor='#FFF9C4',edgecolor='none'))
ax.text(0.2,2.44,'Information ↑',ha='left',va='top',fontsize=5,color='#888')
ax.text(8.8,2.44,'Cost ↑',ha='right',va='top',fontsize=5,color='#888')

# Panel (i): nRF bars
ax.add_patch(mpatches.Rectangle((0.2,1.2),2.8,1.1,facecolor='#FAFAFA',edgecolor='#CCC',linewidth=0.5))
ax.text(1.6,2.2,'(i) Accuracy (nRF ↓)',ha='center',va='center',fontsize=7,fontweight='bold')
for i,(m,v,c) in enumerate(zip(['FT2','Fusang\nmulti-k','IQ-TREE2\nGTR'],[0.084,0.105,0.147],['#43A047','#1E88E5','#E53935'])):
    ax.add_patch(mpatches.Rectangle((0.4+i*0.8,1.4),0.6,v*6,facecolor=c,alpha=0.8))
    ax.text(0.7+i*0.8,1.35,m,ha='center',va='top',fontsize=5)
ax.text(0.2,1.25,'0',ha='left',va='top',fontsize=5)
ax.text(3.0,1.25,'0.15',ha='right',va='top',fontsize=5)

# Panel (ii): scalability
ax.add_patch(mpatches.Rectangle((3.2,1.2),2.8,1.1,facecolor='#FAFAFA',edgecolor='#CCC',linewidth=0.5))
ax.text(4.6,2.2,'(ii) Scalability',ha='center',va='center',fontsize=7,fontweight='bold')
for i,(n,t) in enumerate(zip([200,500,1000,2000,5000,10000],[2.8,18.9,27,55,399,70])):
    x=3.4+i*0.4; y=1.4+np.log10(t+1)*0.5
    ax.plot(x,y,'o',color='#1E88E5',markersize=4)
    if i%2==0: ax.text(x,y+0.05,str(n),ha='center',va='bottom',fontsize=4)
ax.text(3.2,1.25,'n (taxa)',ha='left',va='top',fontsize=5)
ax.text(5.9,1.25,'10K: 70s',ha='right',va='top',fontsize=5,fontweight='bold',color='#1E88E5')

# Panel (iii): indel robustness
ax.add_patch(mpatches.Rectangle((6.2,1.2),2.8,1.1,facecolor='#FAFAFA',edgecolor='#CCC',linewidth=0.5))
ax.text(7.6,2.2,'(iii) Indel Robustness',ha='center',va='center',fontsize=7,fontweight='bold')
irs=[0.005,0.01,0.02,0.03,0.05]
fn=[0.075,0.078,0.080,0.082,0.085]
iq=[0.090,0.110,0.147,0.180,0.220]
for ir,f,i in zip(irs,fn,iq):
    x=6.4+ir*30
    ax.plot(x,1.4+f*6,'o',color='#1E88E5',markersize=3)
    ax.plot(x,1.4+i*6,'o',color='#E53935',markersize=3)
ax.plot([6.4+r*30 for r in irs],[1.4+f*6 for f in fn],'-',color='#1E88E5',linewidth=1)
ax.plot([6.4+r*30 for r in irs],[1.4+i*6 for i in iq],'-',color='#E53935',linewidth=1)
ax.text(6.4,1.3,'Fusang',ha='left',va='top',fontsize=5,color='#1E88E5')
ax.text(6.4,1.22,'IQ-TREE2',ha='left',va='top',fontsize=5,color='#E53935')
ax.text(6.2,1.1,'indel rate →',ha='left',va='bottom',fontsize=4,color='#888')

# Footer
ax.text(5,0.5,'k-mer methods outperform IQ-TREE2 GTR by 1.8× on indel-rich data (p<0.001)',
         ha='center',va='center',fontsize=8,
         bbox=dict(boxstyle='round,pad=0.3',facecolor='#FFF9C4',edgecolor='#F9A825'))

out = r'D:\系统发育树项目\Fusang\Fusang-main'
os.makedirs(out,exist_ok=True)
fig.savefig(os.path.join(out,'GraphicalAbstract.pdf'),dpi=300,bbox_inches='tight')
fig.savefig(os.path.join(out,'GraphicalAbstract.png'),dpi=300,bbox_inches='tight')
print('Graphical Abstract saved')
