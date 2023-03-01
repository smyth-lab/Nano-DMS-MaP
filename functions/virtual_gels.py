import plotly.express as px


def virtual_gel(df, lanewidth, spacing, max_size, marker, molar=False, normalize=False, num_ticks = 10, tickdist=0, verbose=False):
    barcodes = sorted(np.unique(df["barcode"].values))
    
    
    data = {barcode : df[df["barcode"] == barcode]["seq_length"].values for barcode in barcodes}
    width = (1+len(barcodes))*lanewidth
    
    
    gel = np.zeros((max_size, width))
    
    blocksize = max_size/max_size
    if verbose:
        print("width:", width)
        print("blocksize:", blocksize)
        print("lanewidth:", lanewidth)

    #size-weighted
    if molar == False:
        if normalize == True:
            histograms = {barcode: np.histogram(data[barcode], bins=max_size, range=(0, max_size), weights=data[barcode], density = True)[0] for barcode in barcodes}
        else:
            histograms = {barcode: np.histogram(data[barcode], bins=max_size, range=(0, max_size), weights=data[barcode])[0] for barcode in barcodes}
    else:
        if normalize == True:
            histograms = {barcode: np.histogram(data[barcode], bins=max_size, range=(0, max_size), density = True)[0] for barcode in barcodes}
        else:
            histograms = {barcode: np.histogram(data[barcode], bins=max_size, range=(0, max_size))[0] for barcode in barcodes}
        

    label_position = []
    for i, barcode in enumerate(barcodes):
        #print(histograms[barcode].shape)
        #print(int(lanewidth))
        #print(np.tile(histograms[barcode],(lanewidth,1)).shape)
        band_start = int(lanewidth*(i+1))
        band_end = int(lanewidth*(i+2)-np.floor(lanewidth*spacing))
        this_label_position = int(np.floor((band_start+band_end)/2)-1)
        if verbose:
            print(barcode, "band start", band_start, "band end", band_end, "label_position", this_label_position)
        label_position.append(this_label_position)
        #print(band_start, band_end)
        gel[:,band_start:band_end] = np.tile(histograms[barcode],(band_end-band_start,1)).T
    fig = px.imshow(gel, color_continuous_scale='gray_r', origin="lower", aspect="auto")
    
    if tickdist == 0:
        tickdist = np.floor(max_size/num_ticks)
        tickvals = [tickdist*i for i in range(num_ticks)]
        ticktext = [blocksize*i for i in tickvals]
        fig.update_yaxes(tickmode="array", tickvals = tickvals, ticktext = ticktext)
    else:
        fig.update_yaxes(dtick=tickdist)
    
    fig.update_xaxes(tickmode="array", tickvals = label_position, ticktext = barcodes)
    
    """
    if molar == False:
        
    if normalize == True:
        """
    fig.layout.coloraxis.showscale = False
    
    return fig, histograms