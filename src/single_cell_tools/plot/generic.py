def plot_base(ax,
            title="",title_size=30,
            labels=[],labels_size=[30,30],
            xtick_pos=None,ytick_pos=None,tick_sizes=[10,10],
            xtick_labels=None,ytick_labels=None,tick_label_sizes=[20,20],
            grid=False,
            spines=[True,True,True,True],
            axis=None,
            legend=True,
            legend_size=20,
            legend_pos=None,
            legend_title="",
            legend_fontsize=20,
            ):
    """
    Function that simplifies the process of preparing the axis with labels, titles, annotations...

    **Arguments:**

     - **ax**: Axis object
     - **title=None**: Title of the axis
     - **title_size=30**: Font size of the title axis
     - **labels=None**: xlabel, ylabel
     - **labels_size=[30,30]**: fontsize of xlabel and ylabel
     - **xtick_pos=None**: Position of the xticks. None sets them by defaults from the plot.
     - **ytick_pos=None**: Position of the yticks. None sets them by defaults from the plot.
     - **tick_sizes=[10,10]**: Sized of the x and y ticks
     - **xtick_labels=None**: Labels of the xticks
     - **ytick_labels=None**: Labels of the ytikcs
     - **tick_label_sizes=[20,20]**: Sizes of the x and y ticks
     - **grid=False**: If to show the grid.
     - **spines=[True,True,True,True]**: Visibility of the spines [top,left,bottom,right]
     - **axis=None**: If None, axis by default, otherwise the axis limits of the image.
     - **legend**=True: If to print a legend
     - **legend_size**=20: Legend size
     - **legend_pos**=None: Position of the legend. None sets it by default.
     - **legend_title**="": Title of the legend
     - **legend_fontsize**=20: Fontsize of the legend title
    """
    
    if title != None:
        ax.set_title(title,fontsize=title_size)
    
    if labels != None:
        ax.set_xlabel(labels[0],fontsize=labels_size[0])
        ax.set_ylabel(labels[1],fontsize=labels_size[1])
    
    ax.tick_params(axis="x",size=tick_sizes[0],labelsize=tick_label_sizes[0])
    ax.tick_params(axis="y",size=tick_sizes[1],labelsize=tick_label_sizes[1])
    
    if type(xtick_pos) != type(None):
        ax.set_xticks(xtick_pos)
    if type(ytick_pos) != type(None):
        ax.set_yticks(ytick_pos)
    if type(xtick_labels) != type(None):
        ax.set_xticklabels(xtick_labels)
    if type(ytick_labels) != type(None):
        ax.set_yticklabels(ytick_labels)
        
    ax.grid(grid)
    
    ax.spines['top'].set_visible(spines[0])
    ax.spines['right'].set_visible(spines[1])
    ax.spines['bottom'].set_visible(spines[2])
    ax.spines['left'].set_visible(spines[3])

    if axis != None:
        ax.axis(axis)

    if legend:
        ax.legend(fontsize=legend_size,loc=legend_pos,title=legend_title,title_fontsize=legend_fontsize)
    
    return