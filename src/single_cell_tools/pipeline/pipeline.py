class multiple_args:
    
    def __init__(self, *args):
        self.x = [*args]

def demultiplex_step(step):
    
    if type(step) == multiple_args:
        demultiplexed_steps = []
        for option in step.x:
            demultiplexed_steps = np.append(demultiplexed_steps, demultiplex_step(option))
    else:
        demultiplexed_steps = [step]
        for key in step.keys():
            if type(step[key]) == multiple_args:
                demultiplexed_steps_new = []
                for val in step[key].x:
                    for step_copy in demultiplexed_steps:
                        step_copy[key] = val
                        demultiplexed_steps_new.append(step_copy)
                demultiplexed_steps = demultiplexed_steps_new.copy()
                        
    return demultiplexed_steps

def run_pipeline(adata,steps,steps_performed=None):
    
    if steps_performed == None:
        steps_performed = range(len(steps))
        
    for step in steps_performed:
        for demultiplexed in demultiplex_step(steps[step]):
            if "args" in demultiplexed.keys() and "kwargs" in demultiplexed.keys():
                demultiplexed["function"](adata,*demultiplexed["args"],**demultiplexed["kwargs"])
            elif "args" in demultiplexed.keys():
                demultiplexed["function"](adata,*demultiplexed["args"])
            elif "kwargs" in demultiplexed.keys():
                demultiplexed["function"](adata,**demultiplexed["kwargs"])
            else:
                demultiplexed["function"](adata)
                
            if "hps_method" in demultiplexed.keys():
                if "hps_args" in demultiplexed.keys() and "hps_kwargs" in demultiplexed.keys():
                    demultiplexed["hps_method"](adata,*demultiplexed["hps_args"],**demultiplexed["hps_kwargs"])
                elif "hps_args" in demultiplexed.keys():
                    demultiplexed["hps_method"](adata,*demultiplexed["hps_args"])
                elif "hps_kwargs" in demultiplexed.keys():
                    demultiplexed["hps_method"](adata,**demultiplexed["hps_kwargs"])
                else:
                    demultiplexed["hps_method"](adata)
            
    return