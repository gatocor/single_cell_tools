def vline(ax,pos,**kwargs):
    """
    Function to plot a vertical line in the plot.

    **Arguments:**

     - **ax**: Axis object where to plot the line.
     - **pos**: Position in the xaxis where to plot the line.
     - **kwargs**: Keyword arguments to be sent to the function matplotlib.pyplot.vlines.

    **Returns:**
    Adds to axis a line.
    """

    ax.vlines(pos,ax.get_ybound()[0],ax.get_ybound()[1],**kwargs)

    return

def hline(ax,pos,**kwargs):
    """
    Function to plot a horizontal line in the plot.

    **Arguments:**

     - **ax**: Axis object where to plot the line.
     - **pos**: Position in the yaxis where to plot the line.
     - **kwargs**: Keyword arguments to be sent to the function matplotlib.pyplot.vlines.

    **Returns:**
    Adds to axis a line.
    """

    ax.hlines(pos,ax.get_xbound()[0],ax.get_xbound()[1],**kwargs)

    return