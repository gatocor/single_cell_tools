def plot_base(ax,
            title="",title_size=30,
            labels=[],labels_size=[30,30],
            xtick_pos=None,ytick_pos=None,tick_sizes=[10,10],
            xtick_labels=None,ytick_labels=None,tick_label_sizes=[20,20],
            grid=False,
            spines=[True,True,True,True]):
    """Plot base"""
    
    ax.set_title(title,fontsize=title_size)
    
    ax.set_xlabel(labels[0],fontsize=labels_size[0])
    ax.set_ylabel(labels[1],fontsize=labels_size[1])
    
    ax.tick_params(axis="x",size=tick_sizes[0],labelsize=tick_label_sizes[0])
    ax.tick_params(axis="y",size=tick_sizes[1],labelsize=tick_label_sizes[1])
    
    if xtick_pos != None:
        ax.set_xticks(xtick_pos)
    if ytick_pos != None:
        ax.set_yticks(ytick_pos)
    if xtick_labels != None:
        ax.set_xticklabels(xtick_labels)
    if ytick_labels != None:
        ax.set_yticklabels(ytick_labels)
        
    ax.grid(grid)
    
    ax.spines['top'].set_visible(spines[0])
    ax.spines['right'].set_visible(spines[1])
    ax.spines['bottom'].set_visible(spines[2])
    ax.spines['left'].set_visible(spines[3])
    
    return