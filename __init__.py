# Global import of functions
# from .globimport import *

from .ridgeplot import * # ridgeplot
from .tryumap import * # umap from scaling
from .umap_allobs import * # umap for all obs
from .rankobs import * # Rank obs as you want. Useful for filtering
from .describe import * # Describe an anndata
from .assign_sex import * # Assign M, F or MF
from .map_to_dahlin import * # Map to dahlin
from .umap3dPlotting import * # 3D UMAP using plotly
from .splitplot import * # Split an obs into multiple umaps
from .perma_plot import * # Use matplotlib widget and capture images
from .get_from_raw import * # copies .raw.X to X
from .test import *
from .cellcycle_corr import * # Genes correlated to cell cycle
from .gsea import * # GSEA using gseapy
from . import pl as pl
from .umap_grid import * # Plot umap as a grid
