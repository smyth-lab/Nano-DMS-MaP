import plotly.express as px
import numpy as np

# Virtual gels are essentially read length histograms
def plot_virtual_gel(samples, sample_length_dists, max_read_length = 6000, binsize=1, counts=False, mass = False, norm_mass = True, height = 600, width=600, norm_percentile = 100, zmin=0, zmax=1, dtick=500):

    nbins = int(max_read_length / binsize)
    ticks = np.array(np.arange(1, 1+nbins)*binsize)

    matrix = np.zeros((len(samples), int(max_read_length/binsize)))
    per_sample_norm = np.zeros(len(samples))
    
    sample_length_hists = {}
    for row, sample in enumerate(samples):
        read_length_hist = np.histogram(sample_length_dists[sample], bins = nbins, range=[1,1+max_read_length])[0]
        per_sample_norm[row] = np.sum(sample_length_dists[sample])
        sample_length_hists[sample] = read_length_hist
        matrix[row, :] = read_length_hist
        
    if counts:
        if zmax == 1:
            zmax = np.max(matrix)
        fig = px.imshow(matrix.T, color_continuous_scale='gray_r', origin="lower", aspect="auto", x=samples, y=ticks, zmax=zmax)
        fig.update_layout(title="Virtual gel with read counts", width=width, height=height)
        #fig.show()
    
    weights = ticks+int(binsize/2)
    weighted_matrix = weights*matrix
    
    if mass:

        fig = px.imshow(weighted_matrix.T, color_continuous_scale='gray_r', origin="lower", aspect="auto", x=samples, y=ticks)
        fig.update_layout(title="Virtual gel normalized by molecular mass (i.e. weight)", width=width, height=height)
        #fig.show()
    
    if norm_mass:
        weighted_norm_values = np.percentile(weighted_matrix, norm_percentile ,axis=1)
        norm_weighted_matrix = weighted_matrix.T / weighted_norm_values

        fig = px.imshow(norm_weighted_matrix, color_continuous_scale='gray_r', origin="lower", aspect="auto", x=samples, y=ticks, zmin=zmin, zmax=zmax)
        fig.update_layout(title="Virtual gel normalized by molecular mass (i.e. weight) and between samples", width=width, height=height)
        fig.update_yaxes(dtick=dtick)
        #fig.show()

    return fig